// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

//! Read a set of FASTQs, convert into an Iterator over ReadPairs.

use std::path::Path;
use failure::Error;

use fastq::{self, RecordRefIter};
use read_pair::{ReadPair, WhichRead};

use std::io::BufRead;
use utils;


/// A set of corresponding FASTQ representing the different read components from a set of flowcell 'clusters'
#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct InputFastqs {
    pub r1: String,
    pub r2: Option<String>,
    pub i1: Option<String>,
    pub i2: Option<String>,
    pub r1_interleaved: bool,
}

/// Read from a parallel set of FASTQ files.
/// Supports any combination of R1/R2/I1/I2 read files,
/// as well as a interleaved R1/R2 file.
pub struct ReadPairIter {
    iters: [Option<RecordRefIter<Box<BufRead>>>; 4],
    // Each input file can interleave up to 2 -- declare those here
    r1_interleaved: bool,
}

impl ReadPairIter {
    /// Open a `ReadPairIter` given a `InputFastqs` describing a set of FASTQ files
    /// for the available parts of a read.
    pub fn from_fastq_files(input_fastqs: InputFastqs) -> Result<ReadPairIter, Error> {
        Self::new(Some(input_fastqs.r1), input_fastqs.r2, input_fastqs.i1, input_fastqs.i2, input_fastqs.r1_interleaved)
    }

    fn new<P: AsRef<Path>>(r1: Option<P>, r2: Option<P>, i1: Option<P>, i2: Option<P>, r1_interleaved: bool) -> Result<ReadPairIter, Error> {

        let mut iters = [None, None, None, None];

        for (idx, r) in [r1,r2,i1,i2].into_iter().enumerate() {
            match r {
                &Some(ref p) => {
                    let rdr = utils::open_with_gz(p)?;
                    let parser = fastq::Parser::new(rdr);
                    iters[idx] = Some(parser.ref_iter());
                },
                _ => (),
            }
        }

        Ok(ReadPairIter{ iters, r1_interleaved })
    }

    fn get_next(&mut self) -> Result<Option<ReadPair>, Error> {

        let mut iter_ended = false;
        let mut rp = ReadPair::empty();

        for (idx, iter_opt) in self.iters.iter_mut().enumerate() {
            match iter_opt {
                &mut Some(ref mut iter) => {
                    iter.advance()?;
                    {
                        let record = iter.get();
                        if record.is_none() {
                            iter_ended = true;
                        }

                        let which = WhichRead::read_types()[idx];
                        record.map(|r| rp.push_read(&r, which));
                    }

                    // If R1 is interleaved, read another entry
                    // and store it as R2
                    if idx == 0 && self.r1_interleaved {
                        iter.advance()?;
                        let record = iter.get();
                        if record.is_none() {
                            iter_ended = true;
                        }
                        let which = WhichRead::read_types()[idx+1];
                        record.map(|r| rp.push_read(&r, which));
                    }
                },
                &mut None => ()
            }
        }
        // FIXME -- should we error out if the files have different lengths?
        if iter_ended {
            return Ok(None);
        }

        Ok(Some(rp))
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