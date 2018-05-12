// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.


//! Read FASTQs write to shardio buckets by barcode while counting barcodes. 
//! Reads without a correct barcode go to a special 'unbarcoded' bucket.
//! First-pass read subsampling is applied, and bc subsampling is applied.
//! Then go back over the unbarcoded bucket & correct barcodes, using the 
//!! BC counts.  Finally sort each bucket and write out FASTH files.

#![allow(dead_code)]

use std::path::{Path};
use std::collections::HashMap;
use std::hash::{Hash};

use fxhash;
use shardio::{ShardWriter, ShardSender, ShardReaderSet, SortKey};

use rayon::prelude::*;
use utils;

use rand::XorShiftRng;
use rand::Rng;

use serde::ser::Serialize;
use serde::de::DeserializeOwned;

use {Barcode, HasBarcode, FastqProcessor};
use barcode::{BarcodeCorrector};
use std::marker::PhantomData;

use failure::Error;

const NUM_READS_BUFFER: usize = 2048;

struct BcSort;

impl<T> SortKey<T,Barcode> for BcSort where T: HasBarcode {
    fn sort_key(v: &T) -> Barcode {
        v.barcode()
    }
}

fn reduce_counts<K: Hash + Eq>(mut v1: HashMap<K, u32>, mut v2: HashMap<K, u32>) -> HashMap<K, u32> {
    for (k,v) in v2.drain() {
        let counter = v1.entry(k).or_insert(0);
        *counter += v;
    }

    v1
}

fn reduce_counts_err<K: Hash + Eq>(v1: Result<HashMap<K, u32>, Error>, v2: Result<HashMap<K, u32>, Error>) -> Result<HashMap<K, u32>, Error> {
    match (v1, v2) {
        (Ok(m1), Ok(m2)) => Ok(reduce_counts(m1, m2)),
        (Err(e1), _) => Err(e1),
        (_, Err(e2)) => Err(e2), 
    }
}


pub struct SortAndCorrect<ProcType, ReadType> {
    fastq_inputs: Vec<ProcType>,
    phantom: PhantomData<ReadType>,
}

impl<ProcType, ReadType> SortAndCorrect<ProcType, ReadType> where
    ProcType: FastqProcessor<ReadType> + Send + Sync + Clone,
    ReadType: 'static + HasBarcode + Serialize + DeserializeOwned + Send + Sync, 
{
    pub fn new(fastq_inputs: Vec<ProcType>) -> SortAndCorrect<ProcType, ReadType> {
        SortAndCorrect {
            fastq_inputs,
            phantom: PhantomData
        }
    }

    #[inline(never)]
    pub fn sort_bcs(&self,
        read_path: impl AsRef<Path>,
        bc_count_path: impl AsRef<Path>) -> Result<(), Error> {

        let ref writer = ShardWriter::<ReadType, Barcode, BcSort>::new(read_path, 256, 8192, 1<<22);

        // FIXME - configure thread consumption
        let fq_w_senders: Vec<(ProcType, ShardSender<ReadType>)> = 
        self.fastq_inputs.iter().cloned().map(|fq| (fq, writer.get_sender())).collect();

        let bc_counts = fq_w_senders.into_par_iter().map(|(fq, sender)| {
            let desc = fq.description();
            info!("start processing: {}", desc); 

            let counts = self.process_fastq(fq, sender)?;
            info!("done processing: {}", desc);
            Ok(counts)
        }).reduce(|| Ok(HashMap::default()), reduce_counts_err)?;

        utils::write_obj(&bc_counts, bc_count_path.as_ref()).expect("error writing BC counts");
        Ok(())
    }

    fn process_fastq<'a>(&self,
        chunk: ProcType,
        mut bc_sender: ShardSender<ReadType>) -> Result<HashMap<Barcode, u32>, Error> {

        // Counter for BCs
        let mut counts = HashMap::new();

        // Always subsample deterministically
        let mut rand = XorShiftRng::new_unseeded();
        // below: 1.0-r "inverts" the sense of the fraction so that values near
        // 1.0 don't overflow.  Values near zero are a problem.
        let bc_subsample_thresh = chunk.bc_subsample_rate();
        let read_subsample_rate = chunk.read_subsample_rate();

        let fastqs = chunk.fastq_files();
        let iter = ::fastq_read::ReadPairIter::from_fastq_files(fastqs)?;

        for _r in iter {
            let r = _r?;

            // Read subsampling
            if rand.next_f64() < read_subsample_rate {

                match chunk.process_read(r) {
                    None => (),
                    Some(read_pair) => {

                        // count reads on each barcode state
                        let c = counts.entry(read_pair.barcode()).or_insert(0);
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

    pub fn process_unbarcoded<'a>(
        bc_corrected_reads_fn: &str,
        unbarcoded_read_shards: Vec<String>,
        barcode_validator: &BarcodeCorrector<'a>,
        bc_subsample_rate: Option<f64>) -> HashMap<Barcode, u32> {

        let ref writer: ShardWriter<ReadType, Barcode, BcSort> = 
            ShardWriter::new(Path::new(bc_corrected_reads_fn), 16, 1<<12, 1<<21);

        // reader
        let reader = ShardReaderSet::<ReadType, Barcode, BcSort>::open(&unbarcoded_read_shards);
        // below: 1.0-r "inverts" the sense of the fraction so that values near
        // 1.0 don't overflow.  Values near zero are a problem.
        let bc_subsample_thresh = bc_subsample_rate.map_or(0,
                                    |r| {   assert!( r > 0.0, "zero barcode fraction passed in");
                                            ((1.0-r) * u64::max_value() as f64) as usize } );

        let iter = reader.make_chunks(128).iter().cloned().
            map(|x| (x, writer.get_sender())).collect::<Vec<_>>().
            into_par_iter().
            map(|(x, mut read_sender)| {

            // Counter for BCs
            let mut counts = HashMap::new();

            for mut read_pair in reader.iter() {

                // Correct the BC, if possible
                // FIXME
                //barcode_validator.correct_barcode(read_pair);

                // Read subsampling has already occured in process_fastq -- don't repeat here
                // barcode subsampling
                if fxhash::hash(&read_pair.barcode().sequence()) >= bc_subsample_thresh {
                    // count reads on each (gem_group, bc)
                    let c = counts.entry(read_pair.barcode()).or_insert(0);
                    *c += 1;
                    read_sender.send(read_pair)
                }
            }

            counts
        });
        
        iter.reduce(|| HashMap::default(), reduce_counts)
    }
}