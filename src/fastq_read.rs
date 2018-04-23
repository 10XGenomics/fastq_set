use fastq::RecordRefIter;
use new::ReadPair;
use std::io::Read;
use failure::Error;

/// Read from a parallel set of FASTQ files.
/// TODO: support interleaved R1/R2 files.
pub struct ReadPairIter<R: Read> {
    iters: [Option<RecordRefIter<R>>; 4]
}

impl<R: Read> ReadPairIter<R> {
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

        Ok(Some(ReadPair::new2(read_parts)))
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