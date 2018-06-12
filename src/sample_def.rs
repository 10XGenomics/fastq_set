// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

//! Find FASTQs from SampleDefs which get passed to the input of 10x pipelines.
//! Work in progress: need to implement bcl2fastq filename lookup conventions

use std::path::PathBuf;
use InputFastqs;

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

    pub fn get_fastqs(&self) -> Vec<InputFastqs> {
        vec![]
    }
}
