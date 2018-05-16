// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

//! Process 10x FASTQs
//! 
extern crate flate2;
extern crate ordered_float;
extern crate fastq;
extern crate rust_htslib;
extern crate failure;

#[macro_use]
extern crate serde_derive;
extern crate serde;
extern crate bincode;

#[macro_use]
extern crate log;

extern crate rayon;
extern crate shardio;
extern crate fxhash;
extern crate rand;
extern crate serde_json;

pub mod read_pair;
pub mod read_pair_iter;
pub mod sample_def;

pub mod barcode;
pub mod bc_sort;
pub mod utils;

pub mod dna_read;
pub mod rna_read;


#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct InputFastqs {
    r1: String,
    r2: Option<String>,
    i1: Option<String>,
    i2: Option<String>,
    r1_interleaved: bool,
}


#[derive(Serialize, Deserialize, Clone, Copy, PartialOrd, Ord, PartialEq, Eq, Hash)]
struct SSeq {
    length: u8,
    sequence: [u8; 23],
}


impl SSeq {
    pub fn new(seq: &[u8]) -> SSeq {
        assert!(seq.len() <= 23);

        let mut sequence = [0u8; 23];
        sequence[0..seq.len()].copy_from_slice(&seq);

        SSeq {
            length: seq.len() as u8,
            sequence,
        }
    }

    pub fn seq(&self) -> &[u8] {
        &self.sequence[0..self.length as usize]
    }

    pub fn len(&self) -> usize {
        self.length as usize
    }
}

impl AsRef<[u8]> for SSeq {
    fn as_ref(&self) -> &[u8] {
        self.seq()
    }
}

use std::fmt;
impl fmt::Debug for SSeq {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut s = String::new();
        for pos in 0..self.len() {
            s.push(self.sequence[pos] as char);
        }
        write!(f, "{}", s)
    }
}


/// Represent a (possibly-correct 10x barcode sequence)
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
            valid
        }
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
            valid: true
        }
    }

    #[cfg(test)]
    fn test_invalid(sequence: &[u8]) -> Barcode {
        Barcode {
            gem_group: 1,
            sequence: SSeq::new(sequence),
            valid: false
        }
    }

    pub fn valid(&self) -> bool {
        self.valid
    }

    pub fn sequence(&self) -> &[u8] {
        self.sequence.seq()
    }
}

pub trait HasBarcode {
    fn barcode(&self) -> Barcode;
    fn barcode_qual(&self) -> &[u8];
    fn set_barcode(&mut self, barcode: Barcode);
}

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

pub trait HasUmi {
    fn umi(&self) -> Umi;
    fn correct_umi(&mut self, corrected_umi: &[u8]);
}

pub trait ToBam {
    fn to_bam(&self) -> rust_htslib::bam::Record;
}

pub trait AlignableRead {
    fn alignable_reads(&self) -> (Option<&[u8]>, Option<&[u8]>);
}

pub trait FastqProcessor<ReadType> {
    fn process_read(&self, read: read_pair::ReadPair) -> Option<ReadType>;
    fn description(&self) -> String;
    fn fastq_files(&self) -> InputFastqs;

    fn bc_subsample_rate(&self) -> f64;
    fn read_subsample_rate(&self) -> f64;
}