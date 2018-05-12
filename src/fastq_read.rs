use fastq::{self, RecordRefIter};
use raw::{ReadPair, WhichRead};
use FastqFiles;

use failure::Error;
use std::path::Path;
use std::fs::File;

use std::boxed::Box;
use flate2::read::MultiGzDecoder;
use std::io::{Read, BufReader};


/// Read from a parallel set of FASTQ files.
/// TODO: support interleaved R1/R2 files.
pub struct ReadPairIter {
    iters: [Option<RecordRefIter<Box<Read>>>; 4],
    // Each input file can interleave up to 2 -- declare those here
    r1_interleaved: bool,
}

pub fn open_w_gz<P: AsRef<Path>>(p: P) -> Result<Box<Read>, Error> {
    let r = File::open(p.as_ref())?;

    if p.as_ref().extension().unwrap() == "gz" {
        let gz = MultiGzDecoder::new(r);
        let buf_reader = BufReader::with_capacity(32*1024, gz);
        Ok(Box::new(buf_reader))
    } else {
        let buf_reader = BufReader::with_capacity(32*1024, r);
        Ok(Box::new(buf_reader))
    }
}


impl ReadPairIter {
    pub fn from_fastq_files(fastq_files: FastqFiles) -> Result<ReadPairIter, Error> {
        Self::new_gz(Some(fastq_files.r1), fastq_files.r2, fastq_files.i1, fastq_files.i2, fastq_files.r1_interleaved)
    }

    pub fn new_gz<P: AsRef<Path>>(r1: Option<P>, r2: Option<P>, i1: Option<P>, i2: Option<P>, r1_interleaved: bool) -> Result<ReadPairIter, Error> {

        let mut iters = [None, None, None, None];

        let read_types = WhichRead::read_types();

        for (idx, r) in [r1,r2,i1,i2].into_iter().enumerate() {
            match r {
                &Some(ref p) => {
                    let rdr = open_w_gz(p)?;
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

    fn next(&mut self) -> Option<Result<ReadPair, Error>> {
        match self.get_next() {
            Ok(Some(v)) => Some(Ok(v)),
            Ok(None) => None,
            Err(v) => Some(Err(v)),
        }
    }
}