use fastq::{self, RecordRefIter};
use new::{ReadPair, WhichRead};

use failure::Error;
use std::path::Path;
use std::fs::File;

use std::boxed::Box;
use flate2::read::MultiGzDecoder;
use std::io::{Read, BufReader};


/// Read from a parallel set of FASTQ files.
/// TODO: support interleaved R1/R2 files.
pub struct ReadPairIter<R: Read> {
    iters: [Option<RecordRefIter<R>>; 4],
    // Each input file can interleave up to 2 -- declare those here
    targets: [[Option<WhichRead>; 2]; 4]
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


impl<R: Read> ReadPairIter<R> {

    pub fn new_gz<P: AsRef<Path>>(r1: Option<P>, r2: Option<P>, i1: Option<P>, i2: Option<P>) -> Result<ReadPairIter<Box<Read>>, Error> {

        let mut iters = [None, None, None, None];
        let mut targets = [[None,None]; 4];

        let read_types = WhichRead::read_types();

        for (idx, r) in [r1,r2,i1,i2].into_iter().enumerate() {
            match r {
                &Some(ref p) => {
                    let rdr = open_w_gz(p)?;
                    let parser = fastq::Parser::new(rdr);
                    iters[idx] = Some(parser.ref_iter());
                    targets[idx] = [Some(read_types[idx]), None]
                },
                _ => (),
            }
        }

        Ok(ReadPairIter{ iters, targets })
    }

    pub fn new_interleaved_gz<P: AsRef<Path>>(r1_r2_interleaved: P, i1: Option<P>, i2: Option<P>) -> Result<ReadPairIter<Box<Read>>, Error> {

        let mut iters = [None, None, None, None];
        let mut targets = [[None,None]; 4];

        let read_types = WhichRead::read_types();

        let rdr = open_w_gz(r1_r2_interleaved)?;
        let ra_iter = fastq::Parser::new(rdr).ref_iter();
        iters[0] = Some(ra_iter);
        targets[0] = [Some(WhichRead::R1), Some(WhichRead::R2)];

        for (idx, r) in [i1,i2].into_iter().enumerate() {
            match r {
                &Some(ref p) => {
                    let rdr = open_w_gz(p)?;
                    let iter = fastq::Parser::new(rdr).ref_iter();
                    iters[idx] = Some(iter);
                    targets[idx] = [Some(read_types[idx]), None]
                },
                _ => (),
            }
        }

        Ok(ReadPairIter{ iters, targets })
    }


    fn get_next(&mut self) -> Result<Option<ReadPair>, Error> {

        let mut read_parts = [None, None, None, None];
        let mut iter_ended = false;

        for (idx, iter_opt) in self.iters.iter_mut().enumerate() {
            match iter_opt {
                &mut Some(ref mut iter) => {
                    iter.advance()?;
                    let record = iter.get();
                    if record.is_none() {
                        iter_ended = true;
                    }
                    read_parts[idx] = record;
                },
                &mut None => ()
            }
        }
        // FIXME -- should we error out if the files have different lengths?
        if iter_ended {
            if !read_parts.iter().all(|x| x.is_none()) {
                // One of the FASTQ ended before the other -- raise an error
                return Err(format_err!("fastq ended early"));
            }

            return Ok(None);
        }

        Ok(Some(ReadPair::new(read_parts)))
    }
}

impl<R: Read> Iterator for ReadPairIter<R> {
    type Item = Result<ReadPair, Error>;

    fn next(&mut self) -> Option<Result<ReadPair, Error>> {
        match self.get_next() {
            Ok(Some(v)) => Some(Ok(v)),
            Ok(None) => None,
            Err(v) => Some(Err(v)),
        }
    }
}