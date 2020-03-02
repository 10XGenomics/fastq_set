//! Utilities for finding groups of FASTQ files on disk.

pub mod bcl2fastq;
pub mod bcl_processor;

use failure::Error;
use serde::{Deserialize, Serialize};

use crate::read_pair_iter::InputFastqs;

pub use bcl2fastq::Bcl2FastqDef;
pub use bcl_processor::BclProcessorFastqDef;

/// A method to find a set of `InputFastqs` based on
/// some configuration information held by `self`,
/// and some conventions encoded in the implementing
/// type
pub trait FindFastqs {
    fn find_fastqs(&self) -> Result<Vec<InputFastqs>, Error>;
}

/// A pointer to FASTQ data on disk. Can be encoded in the standard Illumina
/// 'bcl2fastq' naming convention, or in the 10x-specific 'BclProcessor'
/// convention. Use the `find_fastqs()` method to find the concrete
/// `InputFastq` files corresponding to a `FastqDef`.
#[derive(Deserialize, Serialize, Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub enum FastqDef {
    Bcl2Fastq(Bcl2FastqDef),
    BclProcessor(BclProcessorFastqDef),
}

impl FastqDef {
    pub fn bcl2fastq(
        fastq_path: String,
        sample_name: String,
        lanes: Option<Vec<usize>>,
    ) -> FastqDef {
        FastqDef::Bcl2Fastq(Bcl2FastqDef {
            fastq_path,
            sample_name,
            lanes,
        })
    }
    pub fn bcl_processor(
        fastq_path: String,
        sample_indices: Vec<String>,
        lanes: Option<Vec<usize>>,
    ) -> FastqDef {
        FastqDef::BclProcessor(BclProcessorFastqDef {
            fastq_path,
            sample_indices,
            lanes,
        })
    }
}

impl FindFastqs for FastqDef {
    fn find_fastqs(&self) -> Result<Vec<InputFastqs>, Error> {
        match self {
            FastqDef::Bcl2Fastq(d) => d.find_fastqs(),
            FastqDef::BclProcessor(d) => d.find_fastqs(),
        }
    }
}
