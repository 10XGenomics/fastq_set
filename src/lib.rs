// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

//! Tools for working with groups FASTQ of files.
//! Major functionality includes:
//! * Find groups FASTQs (R1/R2/I1/I2) following Illumina filename conventions
//! * Parsing flowcell information from Illumina FASTQ headers
//! * High-speed FASTQ I/O (via the `fastq` crate), with careful validation of FASTQ correctness and good error message.
//! * Containers for FASTQ read-pairs (along with index reads), providing access to 'technical' read components like cell barcode and
//! UMI sequences.
//! * Flexible read trimming inspired by `cutadapt`

#![deny(warnings, unused)]
// Allowed clippy lints
#![allow(
    clippy::range_plus_one,
    clippy::naive_bytecount,
    clippy::option_map_unit_fn
)]

pub mod adapter_trimmer;
pub mod array;
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
mod utils;

use crate::read_pair_iter::{AnyReadPairIter, InputFastqs, ReadPairIter};
pub use crate::squality::SQuality;
pub use crate::sseq::SSeq;
use anyhow::Error;
pub use fastq::OwnedRecord;
pub use fastq::Record;
pub use read_pair::WhichRead;
use read_pair_iter::FastqError;
use serde::{Deserialize, Serialize};

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

pub enum ProcessResult<T> {
    Processed(T),
    Unprocessed {
        read: read_pair::ReadPair,
        reason: String,
    },
}

/// A specification for a group of input FASTQ data, and how to interpret
/// the raw sequences as a assay-specific `ReadType` that can provide access to
/// barcodes, UMIs, and track trimmed bases.
pub trait FastqProcessor {
    type ReadType;

    /// Convert a `ReadPair` representing the raw data from a single read-pair
    /// into an assay-specific `ReadType`
    fn process_read(&self, read: read_pair::ReadPair) -> ProcessResult<Self::ReadType>;

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

    fn iter_background_with_storage(
        &self,
        read_ahead: usize,
        storage: read_pair::ReadPairStorage,
    ) -> Result<FastqProcessorIter<'_, Self>, Error>
    where
        Self: Sized,
    {
        FastqProcessorIter::new_background_with_storage(self, read_ahead, storage)
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
    fn make_read_pair_iter(processor: &'a Processor) -> Result<ReadPairIter, Error> {
        let read_pair_iter = ReadPairIter::from_fastq_files(&processor.fastq_files())?
            .illumina_r1_trim_length(processor.illumina_r1_trim_length())
            .illumina_r2_trim_length(processor.illumina_r2_trim_length())
            .subsample_rate(processor.read_subsample_rate());

        Ok(read_pair_iter)
    }

    pub fn new(processor: &'a Processor) -> Result<Self, Error> {
        let iter = Self::make_read_pair_iter(processor)?;
        let read_pair_iter = AnyReadPairIter::Direct(iter);
        Ok(FastqProcessorIter {
            read_pair_iter,
            processor,
        })
    }

    pub fn new_background(processor: &'a Processor, readahead: usize) -> Result<Self, Error> {
        let iter = Self::make_read_pair_iter(processor)?;

        let bg_iter = background_iterator::BackgroundIterator::new(iter, readahead);
        let read_pair_iter = AnyReadPairIter::Background(bg_iter);
        Ok(FastqProcessorIter {
            read_pair_iter,
            processor,
        })
    }

    pub fn new_background_with_storage(
        processor: &'a Processor,
        readahead: usize,
        storage: read_pair::ReadPairStorage,
    ) -> Result<Self, Error> {
        let iter = Self::make_read_pair_iter(processor)?.storage(storage);

        let bg_iter = background_iterator::BackgroundIterator::new(iter, readahead);
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
        let read_pair_iter = Self::make_read_pair_iter(processor)?.storage(storage);

        let read_pair_iter = AnyReadPairIter::Direct(read_pair_iter);
        Ok(FastqProcessorIter {
            read_pair_iter,
            processor,
        })
    }

    pub fn with_seed(processor: &'a Processor, seed: u64) -> Result<Self, Error> {
        let read_pair_iter = Self::make_read_pair_iter(processor)?.seed(seed);

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
        let read_pair_iter = Self::make_read_pair_iter(processor)?
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
    type Item = Result<ProcessResult<<Processor as FastqProcessor>::ReadType>, Error>;

    /// Iterate over ReadType objects.
    fn next(&mut self) -> Option<Self::Item> {
        self.read_pair_iter.next().map(|rr| {
            rr.map(|read| self.processor.process_read(read))
                .map_err(FastqError::into)
        })
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
