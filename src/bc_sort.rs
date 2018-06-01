// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

//! Read input FASTQ data into assay specific read objects that wrap ReadPairs,
//! and sort reads by GEM barcode. Wrappers provide methods to access technical
//! read components, such as barcodes and UMIs. 

use std::path::{Path};

use fxhash;
use shardio::{ShardWriter, ShardSender, ShardReader, Range, SortKey};
use std::marker::PhantomData;
use fxhash::FxHashMap;
use failure::Error;

use rayon::prelude::*;
use utils;

use rand::XorShiftRng;
use rand::Rng;

use serde::ser::Serialize;
use serde::de::DeserializeOwned;

use {Barcode, HasBarcode, FastqProcessor};
use barcode::{BarcodeCorrector, reduce_counts, reduce_counts_err};
use read_pair_iter::ReadPairIter;


pub struct BcSort;

impl<T> SortKey<T,Barcode> for BcSort where T: HasBarcode {
    fn sort_key(v: &T) -> &Barcode {
        v.barcode()
    }
}


pub struct SortByBc<ProcType, ReadType> {
    fastq_inputs: Vec<ProcType>,
    phantom: PhantomData<ReadType>,
}

impl<ProcType, ReadType> SortByBc<ProcType, ReadType> where
    ProcType: FastqProcessor<ReadType> + Send + Sync + Clone,
    ReadType: 'static + HasBarcode + Serialize + DeserializeOwned + Send + Sync, 
{
    pub fn new(fastq_inputs: Vec<ProcType>) -> SortByBc<ProcType, ReadType> {
        SortByBc {
            fastq_inputs,
            phantom: PhantomData
        }
    }

    #[inline(never)]
    pub fn sort_bcs(&self,
        read_path: impl AsRef<Path>,
        bc_count_path: impl AsRef<Path>) -> Result<(u32, u32), Error> {

        let mut writer = ShardWriter::<ReadType, Barcode, BcSort>::new(read_path, 16, 128, 1<<22)?;

        // FIXME - configure thread consumption
        let fq_w_senders: Vec<(ProcType, ShardSender<ReadType>)> = 
        self.fastq_inputs.iter().cloned().map(|fq| (fq, writer.get_sender())).collect();

        let bc_counts = fq_w_senders.into_par_iter().map(|(fq, sender)| {
            let desc = fq.description();
            info!("start processing: {}", desc); 

            let counts = self.process_fastq(fq, sender)?;
            
            info!("done processing: {}", desc);
            Ok(counts)
        }).reduce(|| Ok(FxHashMap::default()), reduce_counts_err)?;

        writer.finish()?;
        let bc: u32 = bc_counts.iter().filter(|(k,_)| k.valid()).map(|v| v.1).sum();
        let no_bc: u32 = bc_counts.iter().filter(|(k,_)| !k.valid()).map(|v| v.1).sum();
        println!("comb bc: {}, no bc: {}", bc, no_bc);

        utils::write_obj(&bc_counts, bc_count_path.as_ref()).expect("error writing BC counts");
        Ok((bc, no_bc))
    }

    fn process_fastq<'a>(&self,
        chunk: ProcType,
        mut bc_sender: ShardSender<ReadType>) -> Result<FxHashMap<Barcode, u32>, Error> {

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
                    Some(read_pair) => {

                        // count reads on each barcode state
                        let c = counts.entry(*read_pair.barcode()).or_insert(0);
                        *c += 1;

                        match read_pair.barcode().valid() {
                            true => {
                                // barcode subsampling
                                let hash = fxhash::hash(&read_pair.barcode().sequence()); 
                                if hash >= bc_subsample_thresh as usize {
                                    bc_sender.send(read_pair)
                                }
                            },
                            false => {
                                // don't bc subsample these reads yet -- that will
                                // happen in the 2nd pass in 'process_unbarcoded'
                                bc_sender.send(read_pair)
                            },
                        }
                    }
                }
            }
        }
        Ok(counts)
    }
}


pub struct CorrectBcs<ReadType> {
    reader: ShardReader<ReadType, Barcode, BcSort>,
    corrector: BarcodeCorrector,
    bc_subsample_rate: Option<f64>,
    phantom: PhantomData<ReadType>,
}

impl<ReadType> CorrectBcs<ReadType> where 
    ReadType: 'static + HasBarcode + Serialize + DeserializeOwned + Send + Sync {

    pub fn new(reader: ShardReader<ReadType, Barcode, BcSort>,
            corrector: BarcodeCorrector,
            bc_subsample_rate: Option<f64>) -> CorrectBcs<ReadType> {

        CorrectBcs {
            reader,
            corrector,
            bc_subsample_rate,
            phantom: PhantomData,
        }
    }

    pub fn process_unbarcoded(&mut self,
        corrected_reads_fn: impl AsRef<Path>,
    ) -> Result<(u32, u32), Error> {

        let mut writer: ShardWriter<ReadType, Barcode, BcSort> = 
            ShardWriter::new(corrected_reads_fn.as_ref(), 16, 1<<12, 1<<21)?;

        // below: 1.0-r "inverts" the sense of the fraction so that values near
        // 1.0 don't overflow.  Values near zero are a problem.
        let bc_subsample_thresh = 0;
        
        /*= self.bc_subsample_rate.map_or(0,
                                    |r| {   assert!( r > 0.0, "zero barcode fraction passed in");
                                            ((1.0-r) * u64::max_value() as f64) as usize } );
*/


        let chunks = self.reader.make_chunks(8, &Range::ends_at(Barcode::first_valid()));
        let chunk_results = 
            chunks.
            into_iter().
            map(|range| (range, writer.get_sender())).collect::<Vec<_>>().
            into_par_iter().
            map(|(range, mut read_sender)| {

                // Counter for BCs
                let mut counts = FxHashMap::default();

                for mut read_pair in self.reader.iter_range(&range) {

                    // Correct the BC, if possible
                    match self.corrector.correct_barcode(&read_pair.barcode(), read_pair.barcode_qual()) {
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
        let w_count = writer.finish()?;

        let bc: u32 = bc_counts.iter().filter(|(k,_)| k.valid()).map(|v| v.1).sum();
        let no_bc: u32 = bc_counts.iter().filter(|(k,_)| !k.valid()).map(|v| v.1).sum();
        println!("post-corr w: {}, bc: {}, no bc: {}", w_count, bc, no_bc);

        Ok((bc, no_bc))
    }
}
