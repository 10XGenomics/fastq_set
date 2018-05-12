extern crate flate2;
extern crate ordered_float;
extern crate fastq;
extern crate rust_htslib;

#[macro_use] 
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

pub mod raw;
pub mod fastq_read;

pub mod barcode;
pub mod bc_sort;
pub mod utils;

pub mod dna_read;
pub mod rna_read;

use std::path::PathBuf;

#[allow(non_camel_case_types)]
#[derive(Serialize, Deserialize, Clone, Copy, PartialOrd, Ord, PartialEq, Eq, Hash)]
enum FastqMode {
    BCL_PROCESSOR,
    ILMN_BCL2FASTQ,
}

/// Pointer to one logical set of FASTQ from a unique (library, gem_group) tuple
#[derive(Serialize, Deserialize, Clone, PartialOrd, Ord, PartialEq, Eq, Hash)]
pub struct SampleDef {
    fastq_mode: FastqMode,
    gem_group: Option<u16>,
    lanes: Option<Vec<usize>>,
    library_type: String,
    read_path: PathBuf,
    sample_indices: Option<Vec<String>>,
    sample_names: Option<Vec<String>>,
}

impl SampleDef {
    pub fn interleaved(&self) -> bool {
        self.fastq_mode == FastqMode::BCL_PROCESSOR
    }



}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct FastqFiles {
    r1: String,
    r2: Option<String>,
    i1: Option<String>,
    i2: Option<String>,
    r1_interleaved: bool,
}

impl SampleDef {
    pub fn get_fastqs(&self) -> Vec<FastqFiles> {
        vec![]
    }
}

/*
Possible conversions:
SampleDef --> FASTQ filenames --> Iterator of Raw Read Pairs

(RawReadPair, ChemistryDefinition) --> RnaReadPair


*/
#[derive(Serialize, Deserialize, Clone, Copy, PartialOrd, Ord, PartialEq, Eq, Hash, Debug)]
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
}

impl AsRef<[u8]> for SSeq {
    fn as_ref(&self) -> &[u8] {
        self.seq()
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
    
    pub fn valid(&self) -> bool {
        self.valid
    }

    pub fn sequence(&self) -> &[u8] {
        self.sequence.seq()
    }
}

pub trait HasBarcode {
    fn barcode(&self) -> Barcode;
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
    fn process_read(&self, read: raw::ReadPair) -> Option<ReadType>;
    fn description(&self) -> String;
    fn fastq_files(&self) -> FastqFiles;

    fn bc_subsample_rate(&self) -> f64;
    fn read_subsample_rate(&self) -> f64;
}

/*
#[cfg(test)]
mod test {
    use super::*;
    use std::collections::HashMap;

    #[test]
    pub fn test_bcl_processor()
    {
        let sd = r#"
[
        {
            "fastq_mode": "BCL_PROCESSOR",
            "gem_group": null,
            "lanes": null,
            "library_type": null,
            "read_path": "test/bcl_processor/",
            "sample_indices": [
                "GGTTTACT"
            ]
        }
]
        "#;

        let mut wl = HashMap::new();
        wl.insert(Vec::from(b"AAAAA".as_ref()), 1);
        wl.insert(Vec::from(b"AAGAC".as_ref()), 2);
        wl.insert(Vec::from(b"ACGAA".as_ref()), 3);
        wl.insert(Vec::from(b"ACGTT".as_ref()), 4);

        let mut counts = HashMap::new();
        counts.insert((0,1), 100);
        counts.insert((0,2), 11);
        counts.insert((0,3), 2);

        let val = BarcodeValidator {
            max_expected_barcode_errors: 1.0,
            bc_confidence_threshold: 0.95,
            whitelist: &wl,
            bc_counts: &counts,
        };

        // Easy
        assert_eq!(val.correct_barcode(0, b"AAAAA", &vec![66,66,66,66,66]), Some(1));

        // Low quality
        assert_eq!(val.correct_barcode(0, b"AAAAA", &vec![34,34,34,66,66]), None);

        // Trivial correction
        assert_eq!(val.correct_barcode(0, b"AAAAT", &vec![66,66,66,66,40]), Some(1));

        // Pseudo-count kills you
        assert_eq!(val.correct_barcode(0, b"ACGAT", &vec![66,66,66,66,66]), None);

        // Quality help you
        assert_eq!(val.correct_barcode(0, b"ACGAT", &vec![66,66,66,66,40]), Some(3));

        // Counts help you
        assert_eq!(val.correct_barcode(0, b"ACAAA", &vec![66,66,66,66,40]), Some(1));
    }
}
*/