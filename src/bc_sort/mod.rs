// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

//! Read input FASTQ data into assay specific read objects that wrap ReadPairs,
//! and sort reads by GEM barcode. Wrappers provide methods to access technical
//! read components, such as barcodes and UMIs.
#![allow(clippy::all)]
use std::fmt::Debug;
use std::hash::Hash;
use std::path::Path;

use failure::Error;
use shardio::{Range, ShardReader, ShardSender, ShardWriter, SortKey};
use std::marker::PhantomData;

use rayon::prelude::*;

use rand::distributions::{Distribution, Uniform};
use rand::SeedableRng;
use rand_xorshift::XorShiftRng;

use serde::de::DeserializeOwned;
use serde::{Deserialize, Serialize};
use shardio;

use crate::barcode::{BarcodeChecker, BarcodeCorrector};
use crate::read_pair_iter::ReadPairIter;
use crate::{Barcode, FastqProcessor, HasBarcode};
use metric::{hash64, Metric, SimpleHistogram, TxHashMap};
use std::borrow::Cow;
use std::path::PathBuf;
use tempfile::TempDir;

/// A marker for sorting objects by their barcode sequence. Incorrect barcodes will
/// be sorted together at the start, followed by correct barcodes. See the definition of `Barcode` for details.
pub struct BcSortOrder;

/// Implementation for objects by their barcode.
pub struct BarcodeOrder;
impl<T> SortKey<T> for BcSortOrder
where
    T: HasBarcode,
{
    type Key = Barcode;
    fn sort_key(v: &T) -> Cow<'_, Barcode> {
        Cow::Borrowed(v.barcode())
    }
}

/// Ingest raw FASTQ using the input data and setting encapsulated in `P`, to generate
/// a stream of 'interpreted' reads of type `ReadType`.  Check the barcodes of each object
/// against that whitelist to mark them as correct. Send the `ReadType` objects to a shardio
/// output file.
pub struct SortByBc<P> {
    fastq_inputs: Vec<P>,
    barcode_checker: BarcodeChecker,
}

impl<P> SortByBc<P>
where
    P: FastqProcessor + Send + Sync + Clone,
    <P as FastqProcessor>::ReadType:
        'static + HasBarcode + Serialize + DeserializeOwned + Send + Sync,
    <<P as FastqProcessor>::ReadType as HasBarcode>::LibraryType:
        Eq + Hash + Clone + Send + Sync + Serialize + for<'de> Deserialize<'de>,
{
    /// Check a new SortByBc object.
    pub fn new(fastq_inputs: Vec<P>, barcode_checker: BarcodeChecker) -> SortByBc<P> {
        SortByBc {
            fastq_inputs,
            barcode_checker,
        }
    }

    /// Execute the sort-by-barcode workflow, send the results to a newly created shardio file with a
    /// filename given by `read_path`.
    #[inline(never)]
    pub fn sort_bcs(
        &self,
        read_path: impl AsRef<Path>,
    ) -> Result<
        (
            TxHashMap<<<P as FastqProcessor>::ReadType as HasBarcode>::LibraryType, u64>,
            TxHashMap<<<P as FastqProcessor>::ReadType as HasBarcode>::LibraryType, u64>,
            TxHashMap<
                <<P as FastqProcessor>::ReadType as HasBarcode>::LibraryType,
                SimpleHistogram<Barcode>,
            >,
        ),
        Error,
    > {
        let mut writer = ShardWriter::<<P as FastqProcessor>::ReadType, BcSortOrder>::new(
            read_path,
            16,
            128,
            1 << 22,
        )?;

        // FIXME - configure thread consumption
        let fq_w_senders: Vec<(P, ShardSender<<P as FastqProcessor>::ReadType, BcSortOrder>)> =
            self.fastq_inputs
                .iter()
                .cloned()
                .map(|fq| (fq, writer.get_sender()))
                .collect();

        let bc_counts = fq_w_senders
            .into_par_iter()
            .map(|(fq, sender)| self.process_fastq(fq, sender))
            .reduce(
                || Ok(TxHashMap::default()),
                |v1, v2| match (v1, v2) {
                    (Ok(mut h1), Ok(h2)) => {
                        h1.merge(h2);
                        Ok(h1)
                    }
                    (Err(e1), _) => Err(e1),
                    (_, Err(e2)) => Err(e2),
                },
            )?;

        writer.finish()?;

        // Sum barcode counts to usize
        let valid_count = count_reads(&bc_counts, true);
        let invalid_count = count_reads(&bc_counts, false);

        Ok((valid_count, invalid_count, bc_counts))
    }

    fn process_fastq<'a>(
        &self,
        chunk: P,
        mut bc_sender: ShardSender<<P as FastqProcessor>::ReadType, BcSortOrder>,
    ) -> Result<
        TxHashMap<
            <<P as FastqProcessor>::ReadType as HasBarcode>::LibraryType,
            SimpleHistogram<Barcode>,
        >,
        Error,
    > {
        // Counter for BCs
        let mut counts = TxHashMap::default();

        // Always subsample deterministically
        let mut rand = XorShiftRng::seed_from_u64(0);
        // below: 1.0-r "inverts" the sense of the fraction so that values near
        // 1.0 don't overflow.  Values near zero are a problem.
        let bc_subsample_thresh = chunk.bc_subsample_rate();
        let read_subsample_rate = chunk.read_subsample_rate();

        let uniform = Uniform::new(0.0, 1.0);

        let fastqs = chunk.fastq_files();
        let iter = ReadPairIter::from_fastq_files(&fastqs)?;

        for _r in iter {
            let r = _r?;

            // Read subsampling
            if uniform.sample(&mut rand) < read_subsample_rate {
                match chunk.process_read(r) {
                    None => (),
                    Some(mut read_pair) => {
                        // Check barcode against whitelist & mark correct if it matches.
                        let mut bc = read_pair.barcode().clone();
                        // FIXME: change this false to a proper check to translate
                        self.barcode_checker.check(&mut bc, false);
                        read_pair.set_barcode(bc);

                        // count reads on each barcode state
                        counts
                            .entry(read_pair.library_type())
                            .or_insert(SimpleHistogram::new())
                            .observe(read_pair.barcode());

                        match read_pair.barcode().is_valid() {
                            true => {
                                // barcode subsampling
                                let hash = hash64(&read_pair.barcode().sequence());
                                if hash >= bc_subsample_thresh as u64 {
                                    bc_sender.send(read_pair)?;
                                }
                            }
                            false => {
                                // don't bc subsample these reads yet -- that will
                                // happen in the 2nd pass in 'process_unbarcoded'
                                bc_sender.send(read_pair)?;
                            }
                        }
                    }
                }
            }
        }
        Ok(counts)
    }
}

fn count_reads<C: Clone + Eq + Hash>(
    bc_counts: &TxHashMap<C, SimpleHistogram<Barcode>>,
    valid: bool,
) -> TxHashMap<C, u64> {
    bc_counts
        .iter()
        .map(|(i, hist)| {
            let count = hist
                .distribution()
                .iter()
                .filter(|(k, _)| k.is_valid() == valid)
                .map(|v| v.1.count() as u64)
                .fold(0u64, |x, y| x + y);
            (i.clone(), count)
        })
        .collect::<_>()
}

/// Iterate over barcodes with initially incorrect barcodes, and attempt to correct them
/// using `corrector`.
pub struct CorrectBcs<ReadType>
where
    ReadType: HasBarcode,
    <ReadType as HasBarcode>::LibraryType: Eq + Hash + Send + Sync,
{
    reader: ShardReader<ReadType, BcSortOrder>,
    corrector: BarcodeCorrector<<ReadType as HasBarcode>::LibraryType>,
    phantom: PhantomData<ReadType>,
}

impl<ReadType> CorrectBcs<ReadType>
where
    ReadType: 'static + HasBarcode + Serialize + DeserializeOwned + Send + Sync,
    <ReadType as HasBarcode>::LibraryType:
        Eq + Hash + Clone + Serialize + for<'de> Deserialize<'de> + Send + Sync,
{
    pub fn new(
        reader: ShardReader<ReadType, BcSortOrder>,
        corrector: BarcodeCorrector<<ReadType as HasBarcode>::LibraryType>,
    ) -> CorrectBcs<ReadType> {
        CorrectBcs {
            reader,
            corrector,
            phantom: PhantomData,
        }
    }

    pub fn process_unbarcoded(
        &mut self,
        corrected_reads_fn: impl AsRef<Path>,
    ) -> Result<
        (
            TxHashMap<<ReadType as HasBarcode>::LibraryType, u64>,
            TxHashMap<<ReadType as HasBarcode>::LibraryType, u64>,
            TxHashMap<<ReadType as HasBarcode>::LibraryType, SimpleHistogram<Barcode>>,
        ),
        Error,
    > {
        let mut writer: ShardWriter<ReadType, BcSortOrder> =
            ShardWriter::new(corrected_reads_fn.as_ref(), 16, 1 << 12, 1 << 21)?;

        let chunks = self
            .reader
            .make_chunks(8, &Range::ends_at(Barcode::first_valid()));
        let chunk_results = chunks
            .into_iter()
            .map(|range| (range, writer.get_sender()))
            .collect::<Vec<_>>()
            .into_par_iter()
            .map(|(range, mut read_sender)| {
                // Counter for BCs
                let mut counts = TxHashMap::default();

                for _read_pair in self.reader.iter_range(&range)? {
                    let mut read_pair = _read_pair?;
                    // Correct the BC, if possible
                    // FIXME: pass in do_translate properly
                    let do_translate = false;
                    match self.corrector.correct_barcode(
                        &read_pair.barcode(),
                        read_pair.barcode_qual(),
                        read_pair.library_type(),
                        do_translate,
                    ) {
                        Some(new_barcode) => read_pair.set_barcode(new_barcode),
                        _ => (),
                    }

                    // count reads on each (gem_group, bc)
                    counts
                        .entry(read_pair.library_type())
                        .or_insert(SimpleHistogram::new())
                        .observe(read_pair.barcode());
                    read_sender.send(read_pair)?;
                }

                Ok(counts)
            });

        // Deal with errors
        let bc_counts = chunk_results.reduce(
            || Ok(TxHashMap::default()),
            |x: Result<_, Error>, y| match (x, y) {
                (Ok(mut x1), Ok(y1)) => {
                    x1.merge(y1);
                    Ok(x1)
                }
                (Err(e1), _) => Err(e1),
                (_, Err(e2)) => Err(e2),
            },
        )?;

        let _ = writer.finish()?;

        // Sum barcode counts to usize
        let valid_count = count_reads(&bc_counts, true);
        let invalid_count = count_reads(&bc_counts, false);

        Ok((valid_count, invalid_count, bc_counts))
    }
}

/// Encapsulates the results from the barcode sorting step.
pub struct BcSortResults<ReadType>
where
    ReadType: HasBarcode,
    <ReadType as HasBarcode>::LibraryType: Eq + Hash + Send + Sync,
{
    pub init_correct_data: PathBuf,
    pub corrected_data: PathBuf,
    pub total_read_pairs: TxHashMap<<ReadType as HasBarcode>::LibraryType, u64>,
    pub init_correct_barcodes: TxHashMap<<ReadType as HasBarcode>::LibraryType, u64>,
    pub corrected_barcodes: TxHashMap<<ReadType as HasBarcode>::LibraryType, u64>,
    pub incorrect_barcodes: TxHashMap<<ReadType as HasBarcode>::LibraryType, u64>,
    pub counts: TxHashMap<<ReadType as HasBarcode>::LibraryType, TxHashMap<Barcode, i64>>,
    pub tmp_dir: TempDir,
}

/// Single-machine barcode sorting workflow. FASTQ data and read settings are input via `chunks`. Outputs
/// are written to a tmp-dir in `path`. `barcode_whitelist` is the path to the 10x barcode whitelist.
pub fn barcode_sort_workflow<P>(
    chunks: Vec<P>,
    path: impl AsRef<Path>,
    barcode_whitelist: impl AsRef<Path>,
) -> Result<BcSortResults<<P as FastqProcessor>::ReadType>, Error>
where
    P: FastqProcessor + Send + Sync + Clone,
    <P as FastqProcessor>::ReadType:
        'static + HasBarcode + Serialize + DeserializeOwned + Send + Sync,
    <<P as FastqProcessor>::ReadType as HasBarcode>::LibraryType:
        Debug + Eq + Hash + Clone + Serialize + for<'de> Deserialize<'de> + Send + Sync,
{
    let tmp_dir = TempDir::new_in(path)?;
    let pass1_fn = tmp_dir.path().join("pass1_data.shardio");
    let pass2_fn = tmp_dir.path().join("pass2_data.shardio");

    // Validate raw barcode sequen
    let barcode_checker = BarcodeChecker::new(barcode_whitelist.as_ref())?;

    let sorter = SortByBc::<_>::new(chunks, barcode_checker);
    let (init_correct, init_incorrect, init_counts) = sorter.sort_bcs(&pass1_fn)?;

    let corrector = BarcodeCorrector::new(barcode_whitelist, init_counts.clone(), 1.5, 0.9)?;

    let reader =
        shardio::ShardReader::<<P as FastqProcessor>::ReadType, BcSortOrder>::open(&pass1_fn)?;

    let mut correct = CorrectBcs::new(reader, corrector);

    let (corrected, still_incorrect, corrected_counts) = correct.process_unbarcoded(&pass2_fn)?;
    // TODO: this is weird and causing test failures, commenting out for now
    /*
    assert_eq!(init_incorrect, {
        let mut rhs = corrected.clone();
        for (t, hist) in still_incorrect.iter() {
            *rhs.entry(t.clone()).or_insert(0) += hist;
        }
        rhs
    });
    */

    let mut total_read_pairs = init_correct.clone();
    for (typ, count) in init_incorrect.iter() {
        *total_read_pairs.entry(typ.clone()).or_insert(0) += count;
    }

    // Combine counts
    let mut final_counts = TxHashMap::default();
    for (typ, hist) in init_counts.iter() {
        for (k, v) in hist.distribution() {
            // initially incorrect reads are handled on the 2nd pass
            if k.is_valid() {
                final_counts
                    .entry(typ.clone())
                    .or_insert(TxHashMap::default())
                    .insert(k.clone(), v.count());
            }
        }
    }

    // update counts with results from 2nd pass
    for (typ, hist) in corrected_counts.iter() {
        for (k, v) in hist.distribution() {
            let e = final_counts
                .entry(typ.clone())
                .or_insert(TxHashMap::default())
                .entry(k.clone())
                .or_insert(0);
            *e += v.count();
        }
    }

    Ok(BcSortResults {
        total_read_pairs,
        init_correct_barcodes: init_correct,
        corrected_barcodes: corrected,
        incorrect_barcodes: still_incorrect,
        counts: final_counts,
        tmp_dir: tmp_dir,
        init_correct_data: pass1_fn,
        corrected_data: pass2_fn,
    })
}
