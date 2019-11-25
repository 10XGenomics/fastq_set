// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

//! Read a set of FASTQs, convert into an Iterator over ReadPairs.

use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};

use crate::read_pair::{MutReadPair, ReadPair, ReadPairStorage, ReadPart, WhichRead};
use fastq::{self, RecordRefIter};

use bytes::{BufMut, BytesMut};
use failure::{format_err, Error, ResultExt};
use rand::distributions::{Distribution, Range};
use rand::{SeedableRng, XorShiftRng};
use std::io::Seek;
use std::io::{BufRead, BufReader, Read};

const GZ_BUF_SIZE: usize = 1 << 16;

/// A set of corresponding FASTQ representing the different read components from a set of flowcell 'clusters'
/// All reads are optional except for R1. For an interleaved R1/R2 file, set the filename in the `r1` field,
/// and set `r1_interleaved = true`.
#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct InputFastqs {
    pub r1: String,
    pub r2: Option<String>,
    pub i1: Option<String>,
    pub i2: Option<String>,
    pub r1_interleaved: bool,
}

impl InputFastqs {
    /// Create a new `InputFastqs` with all paths moved to a new directory
    pub fn change_dir(&self, new_dir: &Path) -> InputFastqs {
        InputFastqs {
            r1: InputFastqs::chng(&self.r1, new_dir),
            r2: self.r2.as_ref().map(|r| InputFastqs::chng(&r, new_dir)),
            i1: self.i1.as_ref().map(|r| InputFastqs::chng(&r, new_dir)),
            i2: self.i2.as_ref().map(|r| InputFastqs::chng(&r, new_dir)),
            r1_interleaved: self.r1_interleaved,
        }
    }

    fn chng(old: &str, new_dir: &Path) -> String {
        let p = Path::new(old);
        let file_name = p.file_name().unwrap();

        let mut new_path = new_dir.to_path_buf();
        new_path.push(file_name);
        new_path.to_string_lossy().to_string()
    }
}

const BUF_SIZE: usize = 4096 * 4;

/// Read sequencing data from a parallel set of FASTQ files.
/// Illumina sequencers typically emit a parallel set of FASTQ files, with one file
/// for each read component taken by the sequencer. Up to 4 reads are possible (R1, R2, I1, and I2).
/// The reader supports any combination of R1/R2/I1/I2 read files,
/// as well as an interleaved R1/R2 file. Supports plain or gzipped FASTQ files, which
/// will be detected based on the filename extension.
pub struct ReadPairIter {
    iters: [Option<RecordRefIter<Box<dyn BufRead>>>; 4],
    paths: [Option<PathBuf>; 4],
    // Each input file can interleave up to 2 -- declare those here
    r1_interleaved: bool,
    buffer: BytesMut,
    rand: XorShiftRng,
    range: Range<f64>,
    subsample_rate: f64,
    storage: ReadPairStorage,
    records_read: [usize; 4],
}

impl ReadPairIter {
    /// Open a `ReadPairIter` given a `InputFastqs` describing a set of FASTQ files
    /// for the available parts of a read.
    pub fn from_fastq_files(input_fastqs: &InputFastqs) -> Result<ReadPairIter, Error> {
        Self::new(
            Some(&input_fastqs.r1),
            input_fastqs.r2.as_ref(),
            input_fastqs.i1.as_ref(),
            input_fastqs.i2.as_ref(),
            input_fastqs.r1_interleaved,
        )
    }

    /// Open a FASTQ file that is uncompressed, gzipped compressed, or lz4 compressed.
    /// The extension of the file is ignored & the filetype is determined by looking
    /// for magic bytes at the of the file
    fn open_fastq(p: impl AsRef<Path>) -> Result<Box<dyn BufRead>, Error> {
        let p = p.as_ref();

        let mut file =
            std::fs::File::open(p).context(format!("Error opening FASTQ file {:?}", p))?;

        let mut buf = vec![0u8; 4];
        file.read_exact(&mut buf[..])?;
        file.seek(std::io::SeekFrom::Start(0))?;

        if &buf[0..2] == &[0x1F, 0x8B] {
            let gz = flate2::read::MultiGzDecoder::new(file);
            let buf_reader = BufReader::with_capacity(GZ_BUF_SIZE, gz);
            Ok(Box::new(buf_reader))
        } else if &buf[0..4] == &[0x04, 0x22, 0x4D, 0x18] {
            let lz = lz4::Decoder::new(file)?;
            let buf_reader = BufReader::with_capacity(GZ_BUF_SIZE, lz);
            Ok(Box::new(buf_reader))
        } else if buf[0] == b'@' {
            let buf_reader = BufReader::with_capacity(32 * 1024, file);
            Ok(Box::new(buf_reader))
        } else {
            let msg = format_err!(
                "FASTQ error: '{:?}' is not appear to be a valid FASTQ. 
            Input FASTQ file must be gzip or lz4 compressed, or must begin with the '@' symbol",
                p
            );
            Err(msg)
        }
    }

    /// Open a (possibly gzip or lz4 compressed) FASTQ file & read some records to confirm the format looks good.
    fn open_fastq_confirm_fmt(p: impl AsRef<Path>) -> Result<Box<dyn BufRead>, Error> {
        let p = p.as_ref();
        let reader = Self::open_fastq(p)?;
        let parser = fastq::Parser::new(reader);

        // make sure we can successfully read some records
        // try and give a useful message if we can't
        let mut iter = parser.ref_iter();

        for _ in 0..10 {
            let read_res = iter.advance();
            read_res.context(format!(
                "Error reading data from FASTQ file '{:?}' - is it corrupted?",
                p
            ))?;

            let rec = iter.get();
            if rec.is_none() {
                break;
            }
        }

        // re-open file so we re-read the initial records
        Self::open_fastq(p)
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

        for (idx, r) in [r1, r2, i1, i2].iter().enumerate() {
            if let Some(ref p) = *r {
                let rdr = Self::open_fastq_confirm_fmt(p)?;
                let parser = fastq::Parser::new(rdr);
                iters[idx] = Some(parser.ref_iter());
                paths[idx] = Some(p.as_ref().to_path_buf());
            }
        }

        let buffer = BytesMut::with_capacity(BUF_SIZE);

        Ok(ReadPairIter {
            paths,
            iters,
            r1_interleaved,
            buffer,
            rand: XorShiftRng::from_seed([42; 16]),
            range: Range::new(0.0, 1.0),
            subsample_rate: 1.0,
            storage: ReadPairStorage::default(),
            records_read: [0; 4],
        })
    }

    pub fn storage(mut self, storage: ReadPairStorage) -> Self {
        self.storage = storage;
        self
    }

    pub fn seed(mut self, seed: [u8; 16]) -> Self {
        self.rand = XorShiftRng::from_seed(seed);
        self
    }

    pub fn subsample_rate(mut self, subsample_rate: f64) -> Self {
        self.subsample_rate = subsample_rate;
        self
    }

    fn get_next(&mut self) -> Result<Option<ReadPair>, Error> {
        // Recycle the buffer if it's almost full.
        if self.buffer.remaining_mut() < 512 {
            self.buffer = BytesMut::with_capacity(BUF_SIZE)
        }

        // need these local reference to avoid borrow checker problem
        let paths = &self.paths;
        let rec_num = &mut self.records_read;

        let mut rp = MutReadPair::empty(&mut self.buffer).storage(self.storage);

        loop {
            let sample = self.range.sample(&mut self.rand) < self.subsample_rate;

            // Track which reader was the first to finish.
            let mut iter_ended = [false; 4];

            for (idx, iter_opt) in self.iters.iter_mut().enumerate() {
                if let Some(ref mut iter) = *iter_opt {
                    iter.advance().with_context(|_| {
                        format!(
                            "FASTQ format error in file {:?} in record starting a line {}",
                            paths[idx].as_ref().unwrap(),
                            rec_num[idx] * 4
                        )
                    })?;

                    {
                        let record = iter.get();
                        // track which reader finished
                        if record.is_none() {
                            iter_ended[idx] = true;
                        }

                        // Check for non-ACGTN characters
                        if let Some(ref rec) = record {
                            if !fastq::Record::validate_dnan(rec) {
                                let msg = format_err!(
                                    "FASTQ contains sequence base with character other than [ACGTN] in file {:?} in record starting on line {}", 
                                    paths[idx].as_ref().unwrap(), rec_num[idx] * 4);
                                return Err(msg);
                            }
                        }

                        if sample {
                            let which = WhichRead::read_types()[idx];
                            record.map(|r| rp.push_read(&r, which));
                        }

                        rec_num[idx] = rec_num[idx] + 1;
                    }

                    // If R1 is interleaved, read another entry
                    // and store it as R2
                    if idx == 0 && self.r1_interleaved && !iter_ended[idx] {
                        iter.advance()?;
                        let record = iter.get();
                        if record.is_none() {
                            // We should only hit this if the FASTQ has an odd
                            // number of records. Throw an error
                            let msg = format_err!("Input FASTQ file {:?} was input as interleaved R1 and R2, but contains an odd number of records .", self.paths[idx].as_ref().unwrap());
                            return Err(msg);
                        }

                        // Check for non-ACGTN characters
                        if let Some(ref rec) = record {
                            if !fastq::Record::validate_dnan(rec) {
                                let msg = format_err!(
                                    "FASTQ contains sequence base with character other than [ACGTN] in file {:?} in record starting at line {}", 
                                    paths[idx].as_ref().unwrap(), rec_num[idx] * 4);
                                return Err(msg);
                            }
                        }

                        if sample {
                            let which = WhichRead::read_types()[idx + 1];
                            record.map(|r| rp.push_read(&r, which));
                        }

                        rec_num[idx] = rec_num[idx] + 1;
                    }
                }
            }

            // check that headers of all reads match
            let mut header_slices = Vec::with_capacity(4);

            let which = [WhichRead::R1, WhichRead::R2, WhichRead::I1, WhichRead::I2];
            for w in 0..4 {
                if let Some(header) = rp.get(which[w], ReadPart::Header) {
                    let prefix = header.split(|x| *x == b' ').next();
                    header_slices.push((w, prefix));
                }
            }

            if header_slices.len() > 0 {
                for i in 1..header_slices.len() {
                    if header_slices[i].1 != header_slices[0].1 {
                        let msg = format_err!("FASTQ header mismatch detected at line {} of input files {:?} and {:?} ",
                                rec_num[0] * 4,
                                self.paths[header_slices[0].0].as_ref().unwrap(),
                                self.paths[header_slices[i].0].as_ref().unwrap()
                            );
                        return Err(msg);
                    }
                }
            }

            // At least one of the readers got to the end -- make sure they all did.
            if iter_ended.iter().any(|x| *x) {
                // Are there any incomplete iterators?
                let any_not_complete = iter_ended
                    .iter()
                    .zip(self.paths.iter())
                    .any(|(ended, path)| path.is_some() && !ended);

                if any_not_complete {
                    // Index of a finished iterator
                    let ended_index = iter_ended.iter().enumerate().find(|(_, v)| **v).unwrap().0;

                    let msg = format_err!(
                        "Input FASTQ file {:?} ended prematurely",
                        self.paths.get(ended_index).unwrap()
                    );
                    return Err(msg);
                } else {
                    return Ok(None);
                }
            }

            if sample {
                return Ok(Some(rp.freeze()));
            }
        }
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
    use file_diff::diff_files;
    use std::fs::File;
    use std::io::Write;

    // Verify that we can parse and write to the identical FASTQ.
    #[test]
    fn test_round_trip() {
        let it = ReadPairIter::new(
            Some("tests/read_pair_iter/good-RA.fastq"),
            None,
            Some("tests/read_pair_iter/good-I1.fastq"),
            Some("tests/read_pair_iter/good-I2.fastq"),
            true,
        )
        .unwrap();

        let _res: Result<Vec<ReadPair>, Error> = it.collect();
        let res = _res.unwrap();

        {
            {
                let mut output = File::create("tests/fastq_round_trip.fastq").unwrap();
                for rec in res {
                    rec.write_fastq(WhichRead::R1, &mut output).unwrap();
                    rec.write_fastq(WhichRead::R2, &mut output).unwrap();
                }
                output.flush().unwrap();
            }

            let mut output = File::open("tests/fastq_round_trip.fastq").unwrap();
            let mut input = File::open("tests/read_pair_iter/good-RA.fastq").unwrap();
            assert!(diff_files(&mut input, &mut output));
        }

        std::fs::remove_file("tests/fastq_round_trip.fastq").unwrap();
    }

    #[test]
    fn test_correct() {
        let it = ReadPairIter::new(
            Some("tests/read_pair_iter/good-RA.fastq"),
            None,
            Some("tests/read_pair_iter/good-I1.fastq"),
            Some("tests/read_pair_iter/good-I2.fastq"),
            true,
        )
        .unwrap();

        let res: Result<Vec<ReadPair>, Error> = it.collect();
        assert!(res.is_ok());
        assert_eq!(res.unwrap().len(), 8);
    }

    #[test]
    fn test_not_gzipped() {
        let it = ReadPairIter::new(
            Some("tests/read_pair_iter/not-gzipped-RA.fastq.gz"),
            None,
            Some("tests/read_pair_iter/good-I1.fastq"),
            Some("tests/read_pair_iter/good-I2.fastq"),
            true,
        )
        .unwrap();

        let res: Result<Vec<ReadPair>, Error> = it.collect();
        assert!(res.is_ok());
    }

    #[test]
    fn test_gzipped() {
        let it = ReadPairIter::new(
            Some("tests/read_pair_iter/good-gzipped-RA.fastq.gz"),
            None,
            Some("tests/read_pair_iter/good-I1.fastq"),
            Some("tests/read_pair_iter/good-I2.fastq"),
            true,
        )
        .unwrap();

        let res: Result<Vec<ReadPair>, Error> = it.collect();
        assert!(res.is_ok());
    }

    #[test]
    fn test_lz4() {
        let it = ReadPairIter::new(
            Some("tests/read_pair_iter/good-lz4-RA.fastq.lz4"),
            None,
            Some("tests/read_pair_iter/good-I1.fastq"),
            Some("tests/read_pair_iter/good-I2.fastq"),
            true,
        )
        .unwrap();

        let res: Result<Vec<ReadPair>, Error> = it.collect();
        assert!(res.is_ok());
    }

    #[test]
    fn test_missing_pair() {
        let it = ReadPairIter::new(
            Some("tests/read_pair_iter/short-RA.fastq"),
            None,
            Some("tests/read_pair_iter/good-I1.fastq"),
            Some("tests/read_pair_iter/good-I2.fastq"),
            true,
        )
        .unwrap();

        let res: Result<Vec<ReadPair>, Error> = it.collect();
        assert!(res.is_err());
    }

    #[test]
    fn test_missing_single_end() {
        let it = ReadPairIter::new(
            Some("tests/read_pair_iter/imbalanced-RA.fastq"),
            None,
            Some("tests/read_pair_iter/good-I1.fastq"),
            Some("tests/read_pair_iter/good-I2.fastq"),
            true,
        )
        .unwrap();

        let res: Result<Vec<ReadPair>, Error> = it.collect();
        assert!(res.is_err());
    }

    #[test]
    fn test_short_i1() {
        let it = ReadPairIter::new(
            Some("tests/read_pair_iter/good-RA.fastq"),
            None,
            Some("tests/read_pair_iter/short-I1.fastq"),
            Some("tests/read_pair_iter/good-I2.fastq"),
            true,
        )
        .unwrap();

        let res: Result<Vec<ReadPair>, Error> = it.collect();
        assert!(res.is_err());
    }

    #[test]
    fn test_bad_char_i1() {
        let it = ReadPairIter::new(
            Some("tests/read_pair_iter/good-RA.fastq"),
            None,
            Some("tests/read_pair_iter/bad-char-I1.fastq"),
            Some("tests/read_pair_iter/good-I2.fastq"),
            true,
        )
        .unwrap();

        let res: Result<Vec<ReadPair>, Error> = it.collect();
        assert!(res.is_err());
    }

    #[test]
    fn test_short_i2() {
        let it = ReadPairIter::new(
            Some("tests/read_pair_iter/good-RA.fastq"),
            None,
            Some("tests/read_pair_iter/good-I1.fastq"),
            Some("tests/read_pair_iter/short-I2.fastq"),
            true,
        )
        .unwrap();

        let res: Result<Vec<ReadPair>, Error> = it.collect();
        assert!(res.is_err());
    }

    #[test]
    fn test_mismatched_header() {
        let it = ReadPairIter::new(
            Some("tests/read_pair_iter/good-RA.fastq"),
            None,
            Some("tests/read_pair_iter/bad-header-I1.fastq"),
            Some("tests/read_pair_iter/good-I2.fastq"),
            true,
        )
        .unwrap();

        let res: Result<Vec<ReadPair>, Error> = it.collect();
        assert!(res.is_err());
    }

    #[test]
    fn test_mismatched_fastq_error() {
        let it = ReadPairIter::new(
            Some("tests/read_pair_iter/good-RA.fastq"),
            None,
            Some("tests/read_pair_iter/bad-format-I1.fastq"),
            Some("tests/read_pair_iter/good-I2.fastq"),
            true,
        );

        // this error gets caught in the opening pre-check
        assert!(it.is_err());
    }

    #[cfg(target_os = "linux")]
    use itertools::Itertools;

    #[cfg(target_os = "linux")]
    fn test_mem_single(every: usize, storage: ReadPairStorage) -> Result<u64, Error> {
        let iter = ReadPairIter::new(
            Some("tests/read_pair_iter/vdj_micro_50k.fastq"),
            None,
            None,
            None,
            true,
        )
        .unwrap()
        .storage(storage);

        let rss_before = psutil::process::Process::new(std::process::id() as i32)?
            .memory()?
            .resident;
        let rp = iter.step_by(every).collect_vec();
        let elements = rp.len();
        let rss_after = psutil::process::Process::new(std::process::id() as i32)?
            .memory()?
            .resident;

        drop(rp);
        let rss_after_drop = psutil::process::Process::new(std::process::id() as i32)?
            .memory()?
            .resident;

        let rss_used = rss_after.saturating_sub(rss_after_drop);
        println!(
            "{:<10} {:<10} {:<12} {:<12} {:<12} {:<12}",
            every,
            elements,
            rss_before / 1024,
            rss_after / 1024,
            rss_used / 1024,
            rss_after_drop / 1024
        );
        Ok(rss_used)
    }

    #[cfg(target_os = "linux")]
    #[test]
    fn test_mem_usage() {
        println!(
            "{:10} {:10} {:12} {:12} {:12} {:12}",
            "Every", "Elements", "Before(kB)", "After(kB)", "Diff(kB)", "Drop(kB)"
        );
        let every = 4;
        let used_rss = test_mem_single(every, ReadPairStorage::PerReadAllocation).unwrap();
        // let used_rss = test_mem_single(every, ReadPairStorage::SharedBuffer).unwrap(); // This fails

        // 50k lines = 6250 reads,
        // With approx 1kB per read the total memory should be ~ 6.25MB / every
        let max_used_rss = 7 * 1024 * 1024 / (every as u64);
        assert!(
            used_rss <= max_used_rss,
            "Used {} bytes whereas max is {}",
            used_rss,
            max_used_rss
        );
    }
}
