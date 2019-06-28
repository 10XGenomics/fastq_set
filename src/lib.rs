// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

//! Containers for FASTQ read-pairs (along with index reads), providing access to 'technical' read components like cell barcode and
//! UMI sequences.
//!

// Allowed clippy lints
#![allow(
    clippy::range_plus_one,
    clippy::naive_bytecount,
    clippy::option_map_unit_fn
)]

#[cfg(test)]
#[macro_use]
extern crate proptest;

#[cfg(test)]
extern crate file_diff;

#[macro_use]
extern crate failure;
extern crate fastq;
extern crate flate2;
extern crate itertools;
extern crate ordered_float;

extern crate glob;
extern crate regex;

#[macro_use]
extern crate serde_derive;
extern crate bincode;
extern crate bytes;
extern crate serde;
extern crate serde_bytes;

extern crate fxhash;
extern crate rand;
extern crate rayon;
extern crate serde_json;
extern crate shardio;
extern crate tempfile;

extern crate bio;
extern crate log;
extern crate lz4;
extern crate metric;

extern crate tenkit;

pub mod read_pair;
pub mod read_pair_iter;
pub mod sample_def;

pub mod barcode;
// bc_sort is an older bc sorting workflow -- it's used in perf tests.
// needs to be reconciled with mod barcode_sort.
pub mod bc_sort;
pub mod barcode_sort;
pub mod metric_utils;
pub mod sseq;
pub mod utils;

pub mod adapter_trimmer;
pub mod dna_read;
pub mod rna_read;

pub use fastq::OwnedRecord;
pub use fastq::Record;

use failure::Error;
use read_pair_iter::{InputFastqs, ReadPairIter};
use sseq::SSeq;

/// Represent a (possibly-corrected) 10x barcode sequence, and it's GEM group
/// FIXME : Should we use the `valid` field for `PartialEq`, `Eq`, `Hash`?
#[derive(Serialize, Deserialize, Clone, Copy, PartialOrd, Ord, PartialEq, Eq, Hash, Debug)]
pub struct Barcode {
    valid: bool,
    gem_group: u16,
    sequence: SSeq,
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

        let mut parts = ss.split("-");
        let bc = parts.next().ok_or_else(|| format_err!("invalid 10x processed barcode: '{}'", ss))?;
        let gg_str = parts.next().ok_or_else(|| format_err!("invalid 10x processed barcode: '{}'", ss))?;

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

    #[cfg(test)]
    fn test_valid(sequence: &[u8]) -> Barcode {
        Barcode {
            gem_group: 1,
            sequence: SSeq::new(sequence),
            valid: true,
        }
    }

    #[cfg(test)]
    fn test_invalid(sequence: &[u8]) -> Barcode {
        Barcode {
            gem_group: 1,
            sequence: SSeq::new(sequence),
            valid: false,
        }
    }

    /// Does this represent a valid whitelist barcode
    pub fn is_valid(&self) -> bool {
        self.valid
    }

    /// Sequence of the (corrected) barcode
    pub fn sequence(&self) -> &[u8] {
        self.sequence.seq()
    }

    /// ASCII string of corrected, GEM group appended form of
    /// barcode, suitable for use in BAM files (CB or BX tags)
    /// For example: "AGCCGATA-1"
    pub fn to_corrected_bytes(&self) -> Vec<u8> {
        let mut v = Vec::with_capacity(self.sequence.len() + 2);
        v.extend(self.sequence());
        v.push(b'-');
        v.extend(format!("{}", self.gem_group).as_bytes());
        v
    }

    pub fn gem_group(&self) -> u16 {
        self.gem_group
    }
}

impl std::fmt::Display for Barcode {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
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
}

/// A trait for objects that carry a UMI sequence.
pub trait HasUmi {
    /// Get a copy of the (possibly corrected) UMI on this read
    fn umi(&self) -> Umi;
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
    /// Convert a `ReadPair` representing the raw data from a single read-pair
    /// into an assay-specific `ReadType`
    fn process_read(&self, read: read_pair::ReadPair) -> Option<Self::ReadType>;

    /// A corresponding set of FASTQ files to read data from.
    fn fastq_files(&self) -> InputFastqs;

    /// Subsampling
    fn bc_subsample_rate(&self) -> f64;
    fn read_subsample_rate(&self) -> f64;

    fn iter(&self) -> Result<FastqProcessorIter<Self>, Error>
    where
        Self: Sized,
    {
        FastqProcessorIter::new(self)
    }

    fn seeded_iter(&self, seed: [u8; 16]) -> Result<FastqProcessorIter<Self>, Error>
    where
        Self: Sized,
    {
        FastqProcessorIter::with_seed(self, seed)
    }

    fn gem_group(&self) -> u16;
}

pub struct FastqProcessorIter<'a, Processor>
where
    Processor: 'a + FastqProcessor,
{
    read_pair_iter: ReadPairIter,
    processor: &'a Processor,
}

impl<'a, Processor> FastqProcessorIter<'a, Processor>
where
    Processor: FastqProcessor,
{
    pub fn new(processor: &'a Processor) -> Result<Self, Error> {
        let read_pair_iter = ReadPairIter::from_fastq_files(processor.fastq_files())?
            .with_subsample_rate(processor.read_subsample_rate());
        Ok(FastqProcessorIter {
            read_pair_iter,
            processor,
        })
    }

    pub fn with_seed(processor: &'a Processor, seed: [u8; 16]) -> Result<Self, Error> {
        let read_pair_iter = ReadPairIter::from_fastq_files(processor.fastq_files())?
            .with_subsample_rate(processor.read_subsample_rate())
            .with_seed(seed);
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

    /// Iterate over ReadType objects
    fn next(&mut self) -> Option<Self::Item> {
        match self.read_pair_iter.next() {
            Some(read_result) => match read_result {
                Ok(read) => Some(Ok(self.processor.process_read(read)?)),
                Err(e) => Some(Err(e)),
            },
            None => None,
        }
    }
}

#[derive(Serialize, Deserialize, Debug, Copy, Clone, PartialEq, Eq)]
pub enum WhichEnd {
    #[serde(rename = "three_prime")]
    ThreePrime,
    #[serde(rename = "five_prime")]
    FivePrime,
}
