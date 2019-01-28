// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

//! Find FASTQs from SampleDefs which get passed to the input of 10x pipelines.
use fxhash::FxHashMap;
use glob::glob;
use regex::Regex;
use std::path::PathBuf;
use tenkit::constants::SAMPLE_INDEX_MAP;
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
        let msg = "Sample indices must be a non-empty list!";
        let indices = self.sample_indices.as_ref().expect(msg);
        if indices.len() == 0 {
            panic!(msg);
        }

        let mut samp_indices: Vec<String> = vec![];
        for index in indices {
            if index == "any" {
                samp_indices.push("*".to_string());
            } else if let Some(seqs) = SAMPLE_INDEX_MAP.get(index.as_str()) {
                samp_indices.extend(seqs.iter().map(|s| s.to_string()));
            } else if is_dna(&index) {
                samp_indices.push(index.clone());
            } else {
                panic!(
                    "Sample index '{}' is not valid. Must be one of: any, SI-<number>, \
                     SI-<plate>-<well-coordinate>, 220<part-number>, or \
                     a nucleotide sequence.",
                    index
                );
            }
        }

        let _types: [&'static str; 3] = ["RA", "I1", "I2"];
        samp_indices
            .iter()
            .flat_map(|si| find_fastqs_10x(&self.read_path, &["RA"], &si, &self.lanes))
            .collect()
    }
}

fn is_dna(seq: &String) -> bool {
    seq.chars().all(|c| "ACGT".contains(c))
}

const MAX_NS: usize = 2;

fn find_fastqs_10x(
    path: &PathBuf,
    read_types: &[&str],
    sample_index: &String,
    lanes: &Option<Vec<usize>>,
) -> Vec<InputFastqs> {
    let si_regex = Regex::new(r".*si-([A-Z]*)_").unwrap();
    let fetch = |pattern: &String| {
        let mut files = vec![];
        for entry in glob(pattern.as_str()).expect("invalid pattern!") {
            match entry {
                Ok(file) => {
                    let f = file
                        .to_str()
                        .expect("non-unicode characters in FASTQ path!");
                    let m = si_regex.captures(f).expect("cannot find sample index!");
                    if m.get(1)
                        .unwrap()
                        .as_str()
                        .chars()
                        .filter(|c| *c == 'N')
                        .count()
                        <= MAX_NS
                    {
                        files.push(file.clone());
                    }
                }
                Err(e) => panic!("unreadable: {:?}", e),
            }
        }
        files.sort_unstable();
        files
    };
    let si_glob = sample_index
        .chars()
        .map(|c| format!("[{}N]", c))
        .collect::<Vec<String>>()
        .join("");
    let prefix = path.join("read-");
    let prefix = prefix
        .to_str()
        .expect("non-unicode characters in FASTQ path!");
    match lanes {
        Some(lanes_) => {
            for lane in lanes_ {
                let mut t2f = FxHashMap::default();
                for read_type in read_types {
                    let pattern = format!(
                        "{}{}_si-{}_lane-{:03}[_\\-]*.fastq*",
                        prefix, read_type, si_glob, lane
                    );
                    let files = fetch(&pattern);
                    t2f.insert(read_type, files);
                }
            }
        }
        None => {
            let mut t2f = FxHashMap::default();
            for read_type in read_types {
                let pattern = format!("{}{}_si-{}*.fastq*", prefix, read_type, si_glob);
                let files = fetch(&pattern);
                t2f.insert(read_type, files);
            }
        }
    }
    vec![]
}

fn find_fastqs_bcl2fastq(
    _path: &PathBuf,
    _read_types: &[&str],
    _sample: &Option<String>,
    _lanes: &Option<Vec<usize>>,
) -> Vec<InputFastqs> {
    vec![]
}
