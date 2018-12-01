// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

//! Read a set of FASTQs, convert into an Iterator over ReadPairs.

use failure::Error;
use std::path::{Path, PathBuf};

use fastq::{self, RecordRefIter};
use read_pair::{MutReadPair, ReadPair, WhichRead};

use bytes::{BytesMut, BufMut};
use failure::format_err;
use std::io::BufRead;
use utils;

/// A set of corresponding FASTQ representing the different read components from a set of flowcell 'clusters'
/// All reads are optional except for R1. For an interleaved R1/R2 file, set the filename in the `r1` field, 
/// and set `r1_interleaved = true`.
#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct InputFastqs {
    pub r1: String,
    pub r2: Option<String>,
    pub i1: Option<String>,
    pub i2: Option<String>,
    pub r1_interleaved: bool,
}

const BUF_SIZE: usize = 4096 * 4;

/// Read sequencing data from a parallel set of FASTQ files.
/// Illumina sequencers typically emit a parallel set of FASTQ files, with one file
/// for each read component taken by the sequencer. Up to 4 reads are possible (R1, R2, I1, and I2).
/// The reader supports any combination of R1/R2/I1/I2 read files, 
/// as well as an interleaved R1/R2 file. Supports plain or gzipped FASTQ files, which
/// will be detected based on the filename extension.
pub struct ReadPairIter {
    iters: [Option<RecordRefIter<Box<BufRead>>>; 4],
    paths: [Option<PathBuf>; 4],
    // Each input file can interleave up to 2 -- declare those here
    r1_interleaved: bool,
    buffer: BytesMut,

}

impl ReadPairIter {
    /// Open a `ReadPairIter` given a `InputFastqs` describing a set of FASTQ files
    /// for the available parts of a read.
    pub fn from_fastq_files(input_fastqs: InputFastqs) -> Result<ReadPairIter, Error> {
        Self::new(
            Some(input_fastqs.r1),
            input_fastqs.r2,
            input_fastqs.i1,
            input_fastqs.i2,
            input_fastqs.r1_interleaved,
        )
    }

    /// Open a `ReadPairIter` given of FASTQ files.
    /// For interleaved R1/R2 files, set `r2 = None`, and set 
    /// `r1_interleaved = true`.
    pub fn new<P: AsRef<Path>>(
        r1: Option<P>,
        r2: Option<P>,
        i1: Option<P>,
        i2: Option<P>,
        r1_interleaved: bool,
    ) -> Result<ReadPairIter, Error> {
        let mut iters = [None, None, None, None];
        let mut paths = [None, None, None, None];

        for (idx, r) in [r1, r2, i1, i2].into_iter().enumerate() {
            match r {
                &Some(ref p) => {
                    let rdr = utils::open_with_gz(p)?;
                    let parser = fastq::Parser::new(rdr);
                    iters[idx] = Some(parser.ref_iter());
                    paths[idx] = Some(p.as_ref().to_path_buf());
                }
                _ => (),
            }
        }

        let buffer = BytesMut::with_capacity(BUF_SIZE);

        Ok(ReadPairIter {
            paths,
            iters,
            r1_interleaved,
            buffer,
        })
    }

    fn get_next(&mut self) -> Result<Option<ReadPair>, Error> {

        // Recycle the buffer if it's almost full.
        if self.buffer.remaining_mut() < 512 {
            self.buffer = BytesMut::with_capacity(BUF_SIZE)
        }

        let mut rp = MutReadPair::empty(&mut self.buffer);

        // Track which reader was the first to finish.
        let mut iter_ended = [false; 4];


        for (idx, iter_opt) in self.iters.iter_mut().enumerate() {
            match iter_opt {
                &mut Some(ref mut iter) => {
                    iter.advance()?;
                    {
                        let record = iter.get();

                        // track which reader finished
                        if record.is_none() {
                            iter_ended[idx] = true;
                        }

                        let which = WhichRead::read_types()[idx];
                        record.map(|r| rp.push_read(&r, which));
                    }

                    // If R1 is interleaved, read another entry
                    // and store it as R2
                    if idx == 0 && self.r1_interleaved && !iter_ended[idx] {
                        iter.advance()?;
                        let record = iter.get();
                        if record.is_none() {
                            // We should only hit this if the FASTQ has an odd
                            // number of records. Throw an error
                            let msg = format_err!("Input FASTQ file {:?} was input as interleaved R1 and R2, but contains an odd number of records .", self.paths.get(idx).unwrap());
                            return Err(msg)
                        }
                        let which = WhichRead::read_types()[idx + 1];
                        record.map(|r| rp.push_read(&r, which));
                    }
                }
                &mut None => (),
            }
        }

        // At least one of the readers got to the end -- make sure they all did.
        if iter_ended.iter().any(|x| *x) {

            // Are there any incomplete iterators?
            let any_not_complete = 
                iter_ended.iter()
                .zip(self.paths.iter())
                .any(|(ended, path)| path.is_some() && !ended);

            if any_not_complete {
                // Index of a finished iterator
                let ended_index = 
                    iter_ended
                    .iter()
                    .enumerate()
                    .filter(|(_, v)| **v).next().unwrap().0;

                let msg = format_err!("Input FASTQ file {:?} ended prematurely", self.paths.get(ended_index).unwrap());
                return Err(msg);
            } else {
                return Ok(None);
            }
        }

        Ok(Some(rp.freeze()))
    }
}

impl Iterator for ReadPairIter {
    type Item = Result<ReadPair, Error>;

    /// Iterate over ReadPair objects
    fn next(&mut self) -> Option<Result<ReadPair, Error>> {
        // Convert Result<Option<_>, Error> to Option<Result<_, Error>>
        match self.get_next() {
            Ok(Some(v)) => Some(Ok(v)),
            Ok(None) => None,
            Err(v) => Some(Err(v)),
        }
    }
}


#[cfg(test)]
mod test_read_pair_iter {
    use super::*;

    #[test]
    fn test_correct() {
        let it = ReadPairIter::new(
            Some("tests/read_pair_iter/good-RA.fastq"),
            None,
            Some("tests/read_pair_iter/good-I1.fastq"),
            Some("tests/read_pair_iter/good-I2.fastq"),
            true
        ).unwrap();

        let res : Result<Vec<ReadPair>, Error> = it.collect();
        assert!(res.is_ok());
        assert_eq!(res.unwrap().len(), 8);
    }

    #[test]
    fn test_missing_pair() {
        let it = ReadPairIter::new(
            Some("tests/read_pair_iter/short-RA.fastq"),
            None,
            Some("tests/read_pair_iter/good-I1.fastq"),
            Some("tests/read_pair_iter/good-I2.fastq"),
            true
        ).unwrap();

        let res : Result<Vec<ReadPair>, Error> = it.collect();
        assert!(res.is_err());
    }

    #[test]
    fn test_missing_single_end() {
        let it = ReadPairIter::new(
            Some("tests/read_pair_iter/imbalanced-RA.fastq"),
            None,
            Some("tests/read_pair_iter/good-I1.fastq"),
            Some("tests/read_pair_iter/good-I2.fastq"),
            true
        ).unwrap();

        let res : Result<Vec<ReadPair>, Error> = it.collect();
        assert!(res.is_err());
    }

    #[test]
    fn test_short_i1() {
        let it = ReadPairIter::new(
            Some("tests/read_pair_iter/good-RA.fastq"),
            None,
            Some("tests/read_pair_iter/short-I1.fastq"),
            Some("tests/read_pair_iter/good-I2.fastq"),
            true
        ).unwrap();

        let res : Result<Vec<ReadPair>, Error> = it.collect();
        assert!(res.is_err());
    }

    #[test]
    fn test_short_i2() {
        let it = ReadPairIter::new(
            Some("tests/read_pair_iter/good-RA.fastq"),
            None,
            Some("tests/read_pair_iter/good-I1.fastq"),
            Some("tests/read_pair_iter/short-I2.fastq"),
            true
        ).unwrap();

        let res : Result<Vec<ReadPair>, Error> = it.collect();
        assert!(res.is_err());
    }
}