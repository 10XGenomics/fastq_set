use failure::Error;
use itertools::Itertools;
use regex;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fmt;
use std::path::{Path, PathBuf};

use super::FindFastqs;
use crate::read_pair_iter::InputFastqs;
use crate::sample_index_map::SAMPLE_INDEX_MAP;

/// A selector for a set of FASTQs emitted by the `BCL_PROCESSOR`
/// demultiplexing pipeline used internally at 10x, and in the
/// (now deprecated) `demux` customer command. `find_fastqs` will
/// select the set of FASTQs in `fastq_path` with the matching
/// sample index and lane values.
#[derive(Deserialize, Serialize, Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub struct BclProcessorFastqDef {
    /// Path to demux / bcl_proccesor FASTQ files
    pub fastq_path: String,
    /// Sample index sequences to include
    pub sample_indices: Vec<String>,
    /// Lanes to include. None indicates that all lanes should be included
    pub lanes: Option<Vec<usize>>,
}

impl FindFastqs for BclProcessorFastqDef {
    fn find_fastqs(&self) -> Result<Vec<InputFastqs>, Error> {
        // get all the files in the directory
        let all_fastq_sets = find_flowcell_fastqs(&self.fastq_path)?;
        let mut res = Vec::new();

        for (bcl_proc, fastqs) in all_fastq_sets {
            // require that the SI match to within one N
            if self
                .sample_indices
                .iter()
                .any(|target| target == "*" || match_seqs_with_n(target, &bcl_proc.si, 1))
            {
                // require that the observed lane is in the allowed list, or it's None
                if self
                    .lanes
                    .as_ref()
                    .map_or(true, |lanes| lanes.contains(&bcl_proc.lane))
                {
                    res.push(fastqs)
                }
            }
        }

        res.sort();
        Ok(res)
    }
}

fn match_seqs_with_n(target: &str, observed: &str, ns_allowed: usize) -> bool {
    let target = target.as_bytes();
    let observed = observed.as_bytes();
    if target.len() != observed.len() {
        return false;
    }

    let mut num_ns = 0;

    for i in 0..target.len() {
        if observed[i] == b'N' {
            num_ns += 1;
        } else if observed[i] != target[i] {
            return false;
        }
    }

    num_ns <= ns_allowed
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct BclProcessorFile {
    pub si: String,
    pub lane: usize,
    pub chunk: usize,
    pub read: String,
    pub path: PathBuf,
}

impl fmt::Display for BclProcessorFile {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self.path)
    }
}

/// Group the FASTQ files according to the known 10x sample index sets
/// declared in `crate::sample_index_map::SAMPLE_INDEX_MAP`.
/// Files are returned in the HashMap where the key is
/// the SI set name, if matched to a known SI set, otherwise it's the SI sequence
/// and the values is a vec of the corresponding input filesets.
pub fn group_samples(
    files: Vec<(BclProcessorFile, InputFastqs)>,
    include_unrec: bool,
) -> HashMap<String, Vec<(BclProcessorFile, InputFastqs)>> {
    let mut si_set_name_map: HashMap<&str, &str> = HashMap::new();

    for (id, vals) in SAMPLE_INDEX_MAP.iter() {
        for si in vals {
            si_set_name_map.insert(si, id);
        }
    }

    let mut by_sid = HashMap::new();

    for (f, input) in files {
        use std::borrow::Borrow;
        let si_str: &str = f.si.borrow();

        match si_set_name_map.get(si_str) {
            Some(set_name) => {
                let entry = by_sid.entry(set_name.to_string()).or_insert(Vec::new());
                entry.push((f, input));
            }
            None => {
                if include_unrec {
                    let entry = by_sid.entry(f.si.clone()).or_insert(Vec::new());
                    entry.push((f, input));
                }
            }
        }
    }

    by_sid
}

pub fn find_flowcell_fastqs(
    path: impl AsRef<Path>,
) -> Result<Vec<(BclProcessorFile, InputFastqs)>, Error> {
    let mut res = Vec::new();

    let mut files = get_demux_files(path)?;
    files.sort();

    for (_, files) in &files
        .into_iter()
        .group_by(|(info, _)| (info.si.clone(), info.lane, info.chunk))
    {
        let my_files: Vec<_> = files.collect();

        let ra = my_files
            .clone()
            .into_iter()
            .find(|(info, _)| info.read == "RA")
            .expect("couldn't find RA read");
        let i1 = my_files
            .clone()
            .into_iter()
            .find(|(info, _)| info.read == "I1");
        let i2 = my_files
            .clone()
            .into_iter()
            .find(|(info, _)| info.read == "I2");

        let fastqs = crate::read_pair_iter::InputFastqs {
            r1: ra.1.to_str().unwrap().to_string(),
            r2: None,
            i1: i1.map(|x| x.1.to_str().unwrap().to_string()),
            i2: i2.map(|x| x.1.to_str().unwrap().to_string()),
            r1_interleaved: true,
        };

        let info = my_files[0].0.clone();
        res.push((info, fastqs));
    }

    Ok(res)
}

fn get_demux_files(path: impl AsRef<Path>) -> Result<Vec<(BclProcessorFile, PathBuf)>, Error> {
    let mut res = Vec::new();
    let dir_files = std::fs::read_dir(path)?;

    for f in dir_files {
        let ent = f?;
        let path = ent.path();

        match try_parse(path) {
            None => (),
            Some(v) => {
                res.push(v);
            }
        }
    }

    Ok(res)
}

fn try_parse(f: PathBuf) -> Option<(BclProcessorFile, PathBuf)> {
    match f.file_name() {
        Some(s) => {
            let r = try_parse_bclprocessor_file(s.to_str().unwrap());
            r.map(|v| (v, f))
        }
        None => None,
    }
}

fn try_parse_bclprocessor_file(filename: &str) -> Option<BclProcessorFile> {
    let re = "^read-([RI][A0-9])_si-([^_]+)_lane-([0-9]+)-chunk-([0-9]+).fastq(.gz)?$";
    let re = regex::Regex::new(re).unwrap();

    match re.captures(filename) {
        None => None,
        Some(caps) => Some(BclProcessorFile {
            path: PathBuf::from(filename),
            read: caps.get(1).unwrap().as_str().to_string(),
            si: caps.get(2).unwrap().as_str().to_string(),
            lane: caps.get(3).unwrap().as_str().parse().unwrap(),
            chunk: caps.get(4).unwrap().as_str().parse().unwrap(),
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bcl_fn() {
        let f1 = "read-RA_si-TTAGCGAT_lane-002-chunk-001.fastq.gz";

        let bcl = try_parse_bclprocessor_file(f1).unwrap();

        let truth = BclProcessorFile {
            path: PathBuf::from(f1),
            read: "RA".to_string(),
            si: "TTAGCGAT".to_string(),
            lane: 2,
            chunk: 1,
        };

        assert_eq!(bcl, truth);
    }

    #[test]
    fn load_dir() -> Result<(), Error> {
        let path = "test/filenames/bcl_processor";
        let fastqs = find_flowcell_fastqs(path)?;
        assert_eq!(fastqs.len(), 44);
        Ok(())
    }

    #[test]
    fn query_all_lanes() -> Result<(), Error> {
        let path = "test/filenames/bcl_processor";

        let query = BclProcessorFastqDef {
            fastq_path: path.to_string(),
            sample_indices: vec!["TCGAATGATC".to_string()],
            lanes: None,
        };

        let fqs = query.find_fastqs()?;
        assert_eq!(fqs.len(), 2);
        Ok(())
    }

    #[test]
    fn query_one_lane() -> Result<(), Error> {
        let path = "test/filenames/bcl_processor";

        let query = BclProcessorFastqDef {
            fastq_path: path.to_string(),
            sample_indices: vec!["TCGAATGATC".to_string()],
            lanes: Some(vec![2]),
        };

        let fqs = query.find_fastqs()?;
        assert_eq!(fqs.len(), 1);
        assert_eq!(
            fqs[0].r1,
            "test/filenames/bcl_processor/read-RA_si-TCGAATGATC_lane-002-chunk-001.fastq.gz"
        );
        Ok(())
    }

    #[test]
    fn test_si_any() {
        let bcl_proc = BclProcessorFastqDef {
            fastq_path: "test/filenames/bcl_processor/".to_string(),
            sample_indices: vec!["*".to_string()],
            lanes: None,
        };
        assert_eq!(bcl_proc.find_fastqs().unwrap().len(), 44);
    }
}
