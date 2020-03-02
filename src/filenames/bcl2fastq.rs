use std::path::{Path, PathBuf};

use failure::Error;
use itertools::Itertools;
use lazy_static::lazy_static;
use regex::Regex;
use serde::{Deserialize, Serialize};

use super::FindFastqs;
use crate::read_pair_iter::InputFastqs;

lazy_static! {
    static ref BCL2FASTQ_REGEX: Regex =
        Regex::new(r"^([\w_-]+)_S(\d+)_L(\d+)_([RI][A12])_(\d+).fastq(.gz)?$").unwrap();
}

/// A pointer to a set of FASTQ files on disk,
/// using the Illumina `bcl2fastq` naming conventions.
/// The `find_fastqs` method will find FASTQ files
/// of the form `heart_1k_v3_S1_L002_R2_001.fastq.gz`
/// with an optional `.gz` suffix.
#[derive(Deserialize, Serialize, Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub struct Bcl2FastqDef {
    /// The path where to the demulitplexed FASTQ files
    pub fastq_path: String,

    /// Sample name used for this sample. Should form a prefix of the filename
    pub sample_name: String,

    /// Lanes to include. None indicates that all lanes should be included
    pub lanes: Option<Vec<usize>>,
}

impl FindFastqs for Bcl2FastqDef {
    fn find_fastqs(&self) -> Result<Vec<InputFastqs>, Error> {
        let all_fastqs = find_flowcell_fastqs(&self.fastq_path)?;

        let mut res = Vec::new();
        for (info, fastqs) in all_fastqs {
            if info.sample == self.sample_name
                && self
                    .lanes
                    .as_ref()
                    .map_or(true, |lanes| lanes.contains(&info.lane))
            // require that the observed lane is in the allowed list, or it's None
            {
                res.push(fastqs);
            }
        }

        res.sort();
        Ok(res)
    }
}

/// A parsed representation of an FASTQ file produced by
/// Illumina's bcl2fastq tool.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct IlmnFastqFile {
    pub sample: String,
    pub s: usize,
    pub lane: usize,
    pub chunk: usize,
    pub read: String,
    pub path: PathBuf,
}

/// Parse the ILMN fastq filename to get the read name, lane, read group,
/// and S field. We expect a filename of the form
/// <path>/<prefix>_S0_L001_R1_001.fastq
impl IlmnFastqFile {
    /// Attempt to parse `path` as an Illumina bcl2fastq-produced
    /// FASTQ file.
    pub fn new(path: impl AsRef<Path>) -> Option<IlmnFastqFile> {
        let filename = path.as_ref().file_name()?.to_str();

        if let Some(f) = filename {
            if let Some(cap) = BCL2FASTQ_REGEX.captures(f) {
                let sample = cap.get(1).unwrap().as_str().to_string();
                let s: usize = cap.get(2).unwrap().as_str().parse().unwrap();
                let lane: usize = cap.get(3).unwrap().as_str().parse().unwrap();
                let read = cap.get(4).unwrap().as_str().to_string();
                let chunk: usize = cap.get(5).unwrap().as_str().parse().unwrap();

                let r = Some(IlmnFastqFile {
                    sample,
                    s,
                    lane,
                    read,
                    chunk,
                    path: path.as_ref().into(),
                });

                return r;
            }
        }

        None
    }
}

/// Find all the bcl2fastq FASTQ files present in `path`.
fn get_bcl2fastq_files(path: impl AsRef<Path>) -> Result<Vec<(IlmnFastqFile, PathBuf)>, Error> {
    let mut res = Vec::new();
    let dir_files = std::fs::read_dir(path)?;

    for f in dir_files {
        let path = f?.path();

        match IlmnFastqFile::new(&path) {
            None => (),
            Some(parsed) => {
                res.push((parsed, path));
            }
        }
    }

    Ok(res)
}

/// Find all the sets of bcl2fastq FASTQ files present in `path`. Corresponding R1/R2/I1/I2 files are grouped
/// together and reported in an `InputFastqs`, along with a representative `IlmnFastqFile` struct.
pub fn find_flowcell_fastqs(
    path: impl AsRef<Path>,
) -> Result<Vec<(IlmnFastqFile, InputFastqs)>, Error> {
    let mut res = Vec::new();

    let mut files = get_bcl2fastq_files(path)?;
    files.sort();

    for (_, files) in &files
        .into_iter()
        .group_by(|(info, _)| (info.sample.clone(), info.s, info.lane, info.chunk))
    {
        let my_files: Vec<_> = files.collect();

        let r1 = my_files
            .clone()
            .into_iter()
            .find(|(info, _)| info.read == "R1");
        let r2 = my_files
            .clone()
            .into_iter()
            .find(|(info, _)| info.read == "R2");
        let i1 = my_files
            .clone()
            .into_iter()
            .find(|(info, _)| info.read == "I1");
        let i2 = my_files
            .clone()
            .into_iter()
            .find(|(info, _)| info.read == "I2");

        // We will tolerate a missing R1 file here -- this
        // could be due to a copying error & unrelated to the
        // sample we're interested in, so we don't want to
        // throw an error just yet.
        if r1.is_none() {
            continue;
        }

        let fastqs = InputFastqs {
            r1: r1.unwrap().1.to_str().unwrap().to_string(),
            r2: r2.map(|x| x.1.to_str().unwrap().to_string()),
            i1: i1.map(|x| x.1.to_str().unwrap().to_string()),
            i2: i2.map(|x| x.1.to_str().unwrap().to_string()),
            r1_interleaved: false,
        };

        let info = my_files[0].0.clone();
        res.push((info, fastqs));
    }

    Ok(res)
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_parse() {
        let filename = "heart_1k_v3_S1_L002_R2_001.fastq.gz";
        let r = IlmnFastqFile::new(filename);

        let expected = IlmnFastqFile {
            path: PathBuf::from(filename.to_string()),
            sample: "heart_1k_v3".to_string(),
            s: 1,
            lane: 2,
            read: "R2".to_string(),
            chunk: 1,
        };

        assert_eq!(r.unwrap(), expected);
    }

    #[test]
    fn test_parse_hyphen() {
        let filename = "heart-1k-v3_S1_L002_R2_001.fastq.gz";
        let r = IlmnFastqFile::new(filename);

        let expected = IlmnFastqFile {
            path: PathBuf::from(filename.to_string()),
            sample: "heart-1k-v3".to_string(),
            s: 1,
            lane: 2,
            read: "R2".to_string(),
            chunk: 1,
        };

        assert_eq!(r.unwrap(), expected);
    }

    #[test]
    fn test_bad() {
        let filename = "heart_1k_v3_S1_LA_R2_001.fastq.gz";
        let r = IlmnFastqFile::new(filename);
        assert!(r.is_none());

        let filename = "heart_1k_v3_S1_L002_XX_001.fastq.gz";
        let r = IlmnFastqFile::new(filename);
        assert!(r.is_none());
    }

    #[test]
    fn query_bcl2fastq() -> Result<(), Error> {
        let path = "test/filenames/bcl2fastq";

        let query = Bcl2FastqDef {
            fastq_path: path.to_string(),
            sample_name: "Infected".to_string(),
            lanes: None,
        };

        let fqs = query.find_fastqs()?;
        assert_eq!(fqs.len(), 2);
        assert_eq!(
            fqs[0].r1,
            "test/filenames/bcl2fastq/Infected_S3_L001_R1_001.fastq"
        );
        assert_eq!(
            fqs[1].r1,
            "test/filenames/bcl2fastq/Infected_S3_L002_R1_001.fastq"
        );
        Ok(())
    }

    #[test]
    fn query_bcl2fastq_lanes() -> Result<(), Error> {
        let path = "test/filenames/bcl2fastq";

        let query = Bcl2FastqDef {
            fastq_path: path.to_string(),
            sample_name: "Infected".to_string(),
            lanes: Some(vec![2]),
        };

        let fqs = query.find_fastqs()?;
        assert_eq!(fqs.len(), 1);
        assert_eq!(
            fqs[0].r1,
            "test/filenames/bcl2fastq/Infected_S3_L002_R1_001.fastq"
        );
        Ok(())
    }
}
