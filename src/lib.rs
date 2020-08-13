// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

//! Tools for working with groups FASTQ of files.
//! Major functionality includes:
//! * Find groups FASTQs (R1/R2/I1/I2) following Illumina filename conventions
//! * Parsing flowcell information from Illumina FASTQ headers
//! * High-speed FASTQ I/O (via the `fastq` crate), with careful validation of FASTQ correctness and good error message.
//! * Containers for FASTQ read-pairs (along with index reads), providing access to 'technical' read components like cell barcode and
//! UMI sequences.
//! * Flexible read trimming inspired by `cutadapt`

#![deny(warnings)]
// Allowed clippy lints
#![allow(
    clippy::range_plus_one,
    clippy::naive_bytecount,
    clippy::option_map_unit_fn
)]

use fastq;

pub mod adapter_trimmer;
pub mod background_iterator;
pub mod filenames;
pub mod illumina_header_info;
pub mod metric_utils;
pub mod read_pair;
pub mod read_pair_iter;
pub mod read_pair_writer;
pub mod sample_index_map;
pub mod squality;
pub mod sseq;
pub mod utils;

pub use crate::squality::SQuality;
pub use crate::sseq::SSeq;

pub use fastq::OwnedRecord;
pub use fastq::Record;

use crate::read_pair_iter::{AnyReadPairIter, InputFastqs, ReadPairIter};
use crate::sseq::{HammingIterOpt, SSeqOneHammingIter};
use failure::{format_err, Error};
use itertools::Itertools;
pub use read_pair::WhichRead;
use serde::{Deserialize, Serialize};

/// Represent a (possibly-corrected) 10x barcode sequence, and it's GEM group.
/// Note the natural sort order groups barcodes by GEM group, then whether they valid, then barcode sequence.
#[derive(Serialize, Deserialize, Clone, Copy, PartialOrd, Ord, PartialEq, Eq, Hash, Debug)]
pub struct Barcode {
    pub gem_group: u16,
    pub valid: bool,
    pub sequence: SSeq,
}

impl Barcode {
    pub fn new(gem_group: u16, sequence: &[u8], valid: bool) -> Barcode {
        Barcode {
            gem_group,
            sequence: SSeq::new(sequence),
            valid,
        }
    }

    pub fn from_sequence(seq: &[u8]) -> Result<Barcode, Error> {
        let ss = std::str::from_utf8(seq)?;

        let mut parts = ss.split('-');
        let bc = parts
            .next()
            .ok_or_else(|| format_err!("invalid 10x processed barcode: '{}'", ss))?;
        let gg_str = parts
            .next()
            .ok_or_else(|| format_err!("invalid 10x processed barcode: '{}'", ss))?;

        use std::str::FromStr;
        let gg = u16::from_str(gg_str)?;
        Ok(Barcode::new(gg, bc.as_bytes(), true))
    }

    /// First possible valid barcode value. Use to query valid/invlaid barcode ranges.
    pub fn first_valid() -> Barcode {
        Barcode {
            valid: true,
            gem_group: 0,
            sequence: SSeq::new(b""),
        }
    }

    /// Create a test barcode, with the given sequence, marked as valid
    pub fn test_valid(sequence: &[u8]) -> Barcode {
        Barcode {
            gem_group: 1,
            sequence: SSeq::new(sequence),
            valid: true,
        }
    }

    /// Create a test barcode, with the given sequence, marked as invalid
    pub fn test_invalid(sequence: &[u8]) -> Barcode {
        Barcode {
            gem_group: 1,
            sequence: SSeq::new(sequence),
            valid: false,
        }
    }

    /// Does this represent a valid whitelist barcode
    pub fn is_valid(self) -> bool {
        self.valid
    }

    /// Sequence of the (corrected) barcode
    pub fn sequence(&self) -> &[u8] {
        self.sequence.seq()
    }

    pub fn sseq(self) -> SSeq {
        self.sequence
    }

    /// ASCII string of corrected, GEM group appended form of
    /// barcode, suitable for use in BAM files (CB or BX tags)
    /// For example: "AGCCGATA-1"
    pub fn to_corrected_bytes(self) -> Vec<u8> {
        let mut v = Vec::with_capacity(self.sequence.len() + 2);
        v.extend(self.sequence());
        v.push(b'-');
        v.extend(format!("{}", self.gem_group).as_bytes());
        v
    }

    pub fn gem_group(self) -> u16 {
        self.gem_group
    }

    pub fn one_hamming_iter(self, opt: HammingIterOpt) -> SSeqOneHammingIter {
        self.sequence.one_hamming_iter(opt)
    }
}

impl std::fmt::Display for Barcode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            String::from_utf8(self.to_corrected_bytes()).unwrap()
        )
    }
}

/// A trait for objects that carry a 10x barcode, allowing for querying the barcode,
/// and correcting the barcode.
pub trait HasBarcode {
    fn barcode(&self) -> &Barcode;
    fn barcode_qual(&self) -> &[u8];
    fn set_barcode(&mut self, barcode: Barcode);
    fn raw_bc_seq(&self) -> &[u8];
    fn raw_bc_qual(&self) -> &[u8];
}

/// A trait for reads that may have a sample index.
pub trait HasSampleIndex {
    fn si_seq(&self) -> Option<&[u8]>;
    fn si_qual(&self) -> Option<&[u8]>;
}

/// A container for a read UMI sequence.
#[derive(Serialize, Deserialize, Clone, Copy, PartialOrd, Ord, PartialEq, Eq, Hash, Debug)]
pub struct Umi {
    sequence: SSeq,
}

impl Umi {
    pub fn new(sequence: &[u8]) -> Umi {
        Umi {
            sequence: SSeq::new(sequence),
        }
    }

    pub fn sequence(&self) -> &[u8] {
        self.sequence.seq()
    }

    pub fn sseq(self) -> SSeq {
        self.sequence
    }

    // No N's and not a homopolymer
    pub fn is_valid(self) -> bool {
        let seq = self.sequence.seq();
        let is_homopolymer = seq.iter().tuple_windows().all(|(a, b)| a == b);
        let has_n = seq.iter().any(|&s| s == b'N' || s == b'n');
        !(is_homopolymer || has_n)
    }

    pub fn one_hamming_iter(self, opt: HammingIterOpt) -> UmiOneHammingIter {
        UmiOneHammingIter {
            inner: self.sequence.one_hamming_iter(opt),
        }
    }
}

pub struct UmiOneHammingIter {
    inner: SSeqOneHammingIter,
}
impl Iterator for UmiOneHammingIter {
    type Item = Umi;

    fn next(&mut self) -> Option<Self::Item> {
        self.inner.next().map(|s| Umi { sequence: s })
    }
}

impl std::fmt::Display for Umi {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        std::fmt::Display::fmt(&self.sequence, f)
    }
}

/// A trait for objects that carry a UMI sequence.
pub trait HasUmi {
    /// Get a copy of the (possibly corrected) UMI on this read
    fn umi(&self) -> Umi;
    fn raw_umi(&self) -> Umi;
    fn correct_umi(&mut self, corrected_umi: &[u8]);
}

/// A trait for objects that carry alignable sequence data.
pub trait AlignableReadPair {
    /// The FASTQ header of the underlying Illumina read
    fn header(&self) -> &[u8];

    /// Alignable sequence bases, representing a trimmed version of the input sequence
    fn alignable_sequence(&self) -> (&[u8], &[u8]);

    /// Quality scores corresponding to the alignable sequences.
    fn alignable_quals(&self) -> (&[u8], &[u8]);
}

/// Specifices what BAM tags should be used to encode the non-alignable
/// parts of the read sequence as BAM tags for BAM to FASTQ conversion
pub trait HasBamTags {
    fn tags(&self) -> Vec<([u8; 2], &[u8])>;
}

/// A specification for a group of input FASTQ data, and how to interpret
/// the raw sequences as a assay-specific `ReadType` that can provide access to
/// barcodes, UMIs, and track trimmed bases.
pub trait FastqProcessor {
    type ReadType;

    /// Determine whether to ignore or propagate
    fn ignore_incomplete_reads(&self) -> bool {
        return true;
    }

    /// Convert a `ReadPair` representing the raw data from a single read-pair
    /// into an assay-specific `ReadType`
    fn process_read(&self, read: read_pair::ReadPair) -> Result<Self::ReadType, Error>;

    /// A corresponding set of FASTQ files to read data from.
    fn fastq_files(&self) -> InputFastqs;

    /// Subsampling
    fn bc_subsample_rate(&self) -> f64;
    fn read_subsample_rate(&self) -> f64;

    /// Read trimming
    fn illumina_r1_trim_length(&self) -> Option<usize>;
    fn illumina_r2_trim_length(&self) -> Option<usize>;

    fn iter(&self) -> Result<FastqProcessorIter<'_, Self>, Error>
    where
        Self: Sized,
    {
        FastqProcessorIter::new(self)
    }

    fn iter_background(&self, read_ahead: usize) -> Result<FastqProcessorIter<'_, Self>, Error>
    where
        Self: Sized,
    {
        FastqProcessorIter::new_background(self, read_ahead)
    }

    fn iter_with_storage(
        &self,
        storage: read_pair::ReadPairStorage,
    ) -> Result<FastqProcessorIter<'_, Self>, Error>
    where
        Self: Sized,
    {
        FastqProcessorIter::with_storage(self, storage)
    }

    fn seeded_iter(&self, seed: u64) -> Result<FastqProcessorIter<'_, Self>, Error>
    where
        Self: Sized,
    {
        FastqProcessorIter::with_seed(self, seed)
    }

    fn seeded_iter_with_storage(
        &self,
        seed: u64,
        storage: read_pair::ReadPairStorage,
    ) -> Result<FastqProcessorIter<'_, Self>, Error>
    where
        Self: Sized,
    {
        FastqProcessorIter::with_seed_and_storage(self, seed, storage)
    }

    fn gem_group(&self) -> u16;
}

pub struct FastqProcessorIter<'a, Processor>
where
    Processor: FastqProcessor,
{
    read_pair_iter: AnyReadPairIter,
    processor: &'a Processor,
}

impl<'a, Processor> FastqProcessorIter<'a, Processor>
where
    Processor: FastqProcessor,
{
    pub fn new(processor: &'a Processor) -> Result<Self, Error> {
        let read_pair_iter = ReadPairIter::from_fastq_files(&processor.fastq_files())?
            .illumina_r1_trim_length(processor.illumina_r1_trim_length())
            .illumina_r2_trim_length(processor.illumina_r2_trim_length())
            .subsample_rate(processor.read_subsample_rate());

        let read_pair_iter = AnyReadPairIter::Direct(read_pair_iter);
        Ok(FastqProcessorIter {
            read_pair_iter,
            processor,
        })
    }

    pub fn new_background(processor: &'a Processor, readahead: usize) -> Result<Self, Error> {
        let read_pair_iter = ReadPairIter::from_fastq_files(&processor.fastq_files())?
            .illumina_r1_trim_length(processor.illumina_r1_trim_length())
            .illumina_r2_trim_length(processor.illumina_r2_trim_length())
            .subsample_rate(processor.read_subsample_rate());

        let bg_iter = background_iterator::BackgroundIterator::new(read_pair_iter, readahead);
        let read_pair_iter = AnyReadPairIter::Background(bg_iter);
        Ok(FastqProcessorIter {
            read_pair_iter,
            processor,
        })
    }

    pub fn with_storage(
        processor: &'a Processor,
        storage: read_pair::ReadPairStorage,
    ) -> Result<Self, Error> {
        let read_pair_iter = ReadPairIter::from_fastq_files(&processor.fastq_files())?
            .subsample_rate(processor.read_subsample_rate())
            .illumina_r1_trim_length(processor.illumina_r1_trim_length())
            .illumina_r2_trim_length(processor.illumina_r2_trim_length())
            .storage(storage);

        let read_pair_iter = AnyReadPairIter::Direct(read_pair_iter);
        Ok(FastqProcessorIter {
            read_pair_iter,
            processor,
        })
    }

    pub fn with_seed(processor: &'a Processor, seed: u64) -> Result<Self, Error> {
        let read_pair_iter = ReadPairIter::from_fastq_files(&processor.fastq_files())?
            .illumina_r1_trim_length(processor.illumina_r1_trim_length())
            .illumina_r2_trim_length(processor.illumina_r2_trim_length())
            .subsample_rate(processor.read_subsample_rate())
            .seed(seed);

        let read_pair_iter = AnyReadPairIter::Direct(read_pair_iter);
        Ok(FastqProcessorIter {
            read_pair_iter,
            processor,
        })
    }

    pub fn with_seed_and_storage(
        processor: &'a Processor,
        seed: u64,
        storage: read_pair::ReadPairStorage,
    ) -> Result<Self, Error> {
        let read_pair_iter = ReadPairIter::from_fastq_files(&processor.fastq_files())?
            .illumina_r1_trim_length(processor.illumina_r1_trim_length())
            .illumina_r2_trim_length(processor.illumina_r2_trim_length())
            .subsample_rate(processor.read_subsample_rate())
            .seed(seed)
            .storage(storage);

        let read_pair_iter = AnyReadPairIter::Direct(read_pair_iter);
        Ok(FastqProcessorIter {
            read_pair_iter,
            processor,
        })
    }
}

impl<'a, Processor> Iterator for FastqProcessorIter<'a, Processor>
where
    Processor: FastqProcessor,
{
    type Item = Result<<Processor as FastqProcessor>::ReadType, Error>;

    /// Iterate over ReadType objects. Reads that aren't parsed succesfully by `process_read`
    /// will be skipped if `FastqProcessor.ignore_incomplete_reads == true`, otherwise
    /// they will be returned as `Some(Err(e))`.
    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let next_read = self.read_pair_iter.next();

            match next_read {
                Some(Ok(read)) => {
                    let r = self.processor.process_read(read);
                    match r {
                        Ok(r) => return Some(Ok(r)),
                        Err(e) => {
                            // if we're not ignoring read parsing error, throw the err.
                            if !self.processor.ignore_incomplete_reads() {
                                return Some(Err(e.into()));
                            }
                            // otherwise we'll try again.
                        }
                    };
                }
                Some(Err(e)) => return Some(Err(e.into())),
                None => return None,
            }
        }
    }
}

/// Which end of a transcript reads come from
#[derive(Serialize, Deserialize, Debug, Copy, Clone, PartialEq, Eq)]
pub enum WhichEnd {
    #[serde(rename = "three_prime")]
    ThreePrime,
    #[serde(rename = "five_prime")]
    FivePrime,
}
