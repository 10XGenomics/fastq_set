// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

//! Read input FASTQ data into assay specific read objects that wrap ReadPairs,
//! and sort reads by GEM barcode. Wrappers provide methods to access technical
//! read components, such as barcodes and UMIs.

use std::path::Path;
use std::borrow::Cow;

use failure::Error;
use fxhash;
use fxhash::FxHashMap;
use shardio::{Range, ShardReader, ShardSender, ShardWriter, SortKey};
use std::marker::PhantomData;

use rayon::prelude::*;

use rand::Rng;
use rand::XorShiftRng;

use serde::de::DeserializeOwned;
use serde::ser::Serialize;
use shardio;

use barcode::{reduce_counts, reduce_counts_err, BarcodeChecker, BarcodeCorrector};
use read_pair_iter::ReadPairIter;
use std::path::PathBuf;
use tempfile::TempDir;
use {Barcode, FastqProcessor, HasBarcode};

/// A marker for sorting objects by their barcode sequence. Incorrect barcodes will
/// be sorted together at the start, followed by correct barcodes. See the definition of `Barcode` for details.
pub struct BcSortOrder;

/// Implementation for objects by their barcode.
impl<T> SortKey<T> for BcSortOrder
where
    T: HasBarcode,
{
    type Key = Barcode;
    fn sort_key(v: &T) -> Cow<Barcode> {
        Cow::from(v.barcode())
    }
}

/// Ingest raw FASTQ using the input data and setting encapsulated in `ProcType`, to generate
/// a stream of 'interpreted' reads of type `ReadType`.  Check the barcodes of each object
/// against that whitelist to mark them as correct. Send the `ReadType` objects to a shardio
/// output file.
pub struct SortByBc<ProcType> {
    fastq_inputs: Vec<ProcType>,
    barcode_checker: BarcodeChecker,
}

impl<ProcType> SortByBc<ProcType>
where
    ProcType: FastqProcessor + Send + Sync + Clone,
    <ProcType as FastqProcessor>::ReadType: 'static + HasBarcode + Serialize + DeserializeOwned + Send + Sync,
{
    /// Check a new SortByBc object.
    pub fn new(
        fastq_inputs: Vec<ProcType>,
        barcode_checker: BarcodeChecker,
    ) -> SortByBc<ProcType> {
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
    ) -> Result<(usize, usize, FxHashMap<Barcode, i64>), Error> {
        let mut writer =
            ShardWriter::<<ProcType as FastqProcessor>::ReadType, Barcode, BcSortOrder>::new(read_path, 16, 128, 1 << 22)?;

        // FIXME - configure thread consumption
        let fq_w_senders: Vec<(ProcType, ShardSender<<ProcType as FastqProcessor>::ReadType>)> = self.fastq_inputs
            .iter()
            .cloned()
            .map(|fq| (fq, writer.get_sender()))
            .collect();

        let bc_counts = fq_w_senders
            .into_par_iter()
            .map(|(fq, sender)| {
                info!("start processing: {:?}", fq.fastq_files());
                self.process_fastq(fq, sender)
            })
            .reduce(|| Ok(FxHashMap::default()), reduce_counts_err)?;

        writer.finish()?;

        // Sum barcode counts to usize
        let valid_count = count_reads(&bc_counts, true);
        let invalid_count = count_reads(&bc_counts, false);

        Ok((valid_count, invalid_count, bc_counts))
    }

    fn process_fastq<'a>(
        &self,
        chunk: ProcType,
        mut bc_sender: ShardSender<<ProcType as FastqProcessor>::ReadType>,
    ) -> Result<FxHashMap<Barcode, u32>, Error> {
        // Counter for BCs
        let mut counts = FxHashMap::default();

        // Always subsample deterministically
        let mut rand = XorShiftRng::new_unseeded();
        // below: 1.0-r "inverts" the sense of the fraction so that values near
        // 1.0 don't overflow.  Values near zero are a problem.
        let bc_subsample_thresh = chunk.bc_subsample_rate();
        let read_subsample_rate = chunk.read_subsample_rate();

        let fastqs = chunk.fastq_files();
        let iter = ReadPairIter::from_fastq_files(fastqs)?;

        for _r in iter {
            let r = _r?;

            // Read subsampling
            if rand.gen_range(0.0, 1.0) < read_subsample_rate {
                match chunk.process_read(r) {
                    None => (),
                    Some(mut read_pair) => {
                        // Check barcode against whitelist & mark correct if it matches.
                        let mut bc = read_pair.barcode().clone();
                        self.barcode_checker.check(&mut bc);
                        read_pair.set_barcode(bc);

                        // count reads on each barcode state
                        let c = counts.entry(*read_pair.barcode()).or_insert(0);
                        *c += 1;

                        match read_pair.barcode().is_valid() {
                            true => {
                                // barcode subsampling
                                let hash = fxhash::hash(&read_pair.barcode().sequence());
                                if hash >= bc_subsample_thresh as usize {
                                    bc_sender.send(read_pair)
                                }
                            }
                            false => {
                                // don't bc subsample these reads yet -- that will
                                // happen in the 2nd pass in 'process_unbarcoded'
                                bc_sender.send(read_pair)
                            }
                        }
                    }
                }
            }
        }
        Ok(counts)
    }
}

fn count_reads(bc_counts: &FxHashMap<Barcode, i64>, valid: bool) -> usize {
    bc_counts
        .iter()
        .filter(|(k, _)| k.is_valid() == valid)
        .map(|v| v.1)
        .fold(0usize, |x, y| x + (*y as usize))
}

/// Iterate over barcodes with initially incorrect barcodes, and attempt to correct them
/// using `corrector`.
pub struct CorrectBcs<ReadType: HasBarcode> {
    reader: ShardReader<ReadType, BcSortOrder>,
    corrector: BarcodeCorrector,
    bc_subsample_rate: Option<f64>,
    phantom: PhantomData<ReadType>,
}

impl<ReadType> CorrectBcs<ReadType>
where
    ReadType: 'static + HasBarcode + Serialize + DeserializeOwned + Send + Sync,
{
    pub fn new(
        reader: ShardReader<ReadType, BcSortOrder>,
        corrector: BarcodeCorrector,
        bc_subsample_rate: Option<f64>,
    ) -> CorrectBcs<ReadType> {
        CorrectBcs {
            reader,
            corrector,
            bc_subsample_rate,
            phantom: PhantomData,
        }
    }

    pub fn process_unbarcoded(
        &mut self,
        corrected_reads_fn: impl AsRef<Path>,
    ) -> Result<(usize, usize, FxHashMap<Barcode, i64>), Error> {
        let mut writer: ShardWriter<ReadType, Barcode, BcSortOrder> =
            ShardWriter::new(corrected_reads_fn.as_ref(), 16, 1 << 12, 1 << 21)?;

        // below: 1.0-r "inverts" the sense of the fraction so that values near
        // 1.0 don't overflow.  Values near zero are a problem.
        let bc_subsample_thresh = 0;

        /*= self.bc_subsample_rate.map_or(0,
                                    |r| {   assert!( r > 0.0, "zero barcode fraction passed in");
                                            ((1.0-r) * u64::max_value() as f64) as usize } );
*/

        let chunks = self.reader
            .make_chunks(8, &Range::ends_at(Barcode::first_valid()));
        let chunk_results = chunks
            .into_iter()
            .map(|range| (range, writer.get_sender()))
            .collect::<Vec<_>>()
            .into_par_iter()
            .map(|(range, mut read_sender)| {
                // Counter for BCs
                let mut counts = FxHashMap::default();

                for mut read_pair in self.reader.iter_range(&range) {
                    // Correct the BC, if possible
                    match self.corrector
                        .correct_barcode(&read_pair.barcode(), read_pair.barcode_qual())
                    {
                        Some(new_barcode) => read_pair.set_barcode(new_barcode),
                        _ => (),
                    }

                    // Read subsampling has already occured in process_fastq -- don't repeat here
                    // barcode subsampling
                    if fxhash::hash(&read_pair.barcode().sequence()) >= bc_subsample_thresh {
                        // count reads on each (gem_group, bc)
                        let c = counts.entry(*read_pair.barcode()).or_insert(0);
                        *c += 1;
                        read_sender.send(read_pair)
                    }
                }

                counts
            });

        let bc_counts = chunk_results.reduce(|| FxHashMap::default(), reduce_counts);
        let _ = writer.finish()?;

        // Sum barcode counts to usize
        let valid_count = count_reads(&bc_counts, true);
        let invalid_count = count_reads(&bc_counts, false);

        Ok((valid_count, invalid_count, bc_counts))
    }
}

/// Encapsulates the results from the barcode sorting step.
pub struct BcSortResults {
    init_correct_data: PathBuf,
    corrected_data: PathBuf,
    total_read_pairs: usize,
    init_correct_barcodes: usize,
    corrected_barcodes: usize,
    incorrect_barcodes: usize,
    counts: FxHashMap<Barcode, i64>,
    tmp_dir: TempDir,
}

/// Single-machine barcode sorting workflow. FASTQ data and read settings are input via `chunks`. Outputs
/// are written to a tmp-dir in `path`. `barcode_whitelist` is the path to the 10x barcode whitelist.
pub fn barcode_sort_workflow<P>(
    chunks: Vec<P>,
    path: impl AsRef<Path>,
    barcode_whitelist: impl AsRef<Path>,
) -> Result<BcSortResults, Error>
where
    P: FastqProcessor + Send + Sync + Clone,
    <P as FastqProcessor>::ReadType: 'static + HasBarcode + Serialize + DeserializeOwned + Send + Sync,
{
    let tmp_dir = TempDir::new_in(path)?;
    let pass1_fn = tmp_dir.path().join("pass1_data.shardio");
    let pass2_fn = tmp_dir.path().join("pass2_data.shardio");

    // Validate raw barcode sequen
    let barcode_checker = BarcodeChecker::new(barcode_whitelist.as_ref())?;

    let sorter = SortByBc::<_>::new(chunks, barcode_checker);
    let (init_correct, init_incorrect, init_counts) = sorter.sort_bcs(&pass1_fn)?;

    let corrector = BarcodeCorrector::new(barcode_whitelist, init_counts.clone(), 1.5, 0.9)?;

    let reader = shardio::ShardReader::<<P as FastqProcessor>::ReadType, Barcode, BcSortOrder>::open(&pass1_fn);

    let mut correct = CorrectBcs::new(reader, corrector, Some(1.0f64));

    let (corrected, still_incorrect, corrected_counts) = correct.process_unbarcoded(&pass2_fn)?;
    assert_eq!(init_incorrect, corrected + still_incorrect);

    // Combine counts
    let mut final_counts = FxHashMap::default();
    for (k, v) in init_counts {
        // initially incorrect reads are handled on the 2nd pass
        if k.is_valid() {
            final_counts.insert(k, v);
        }
    }

    // update counts with results from 2nd pass
    for (k, v) in corrected_counts {
        let e = final_counts.entry(k).or_insert(0);
        *e += v;
    }

    Ok(BcSortResults {
        total_read_pairs: init_correct + init_incorrect,
        init_correct_barcodes: init_correct,
        corrected_barcodes: corrected,
        incorrect_barcodes: still_incorrect,
        counts: final_counts,
        tmp_dir: tmp_dir,
        init_correct_data: pass1_fn,
        corrected_data: pass2_fn,
    })
}
