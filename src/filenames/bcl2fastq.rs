use super::FindFastqs;
use crate::filenames::LaneMode;
use crate::filenames::LaneSpec;
use crate::read_pair_iter::InputFastqs;
use failure::Error;
use itertools::Itertools;
use lazy_static::lazy_static;
use metric::TxHashSet;
use regex::Regex;
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};

lazy_static! {
    static ref BCL2FASTQ_REGEX: Regex =
        Regex::new(r"^([\w_-]+)_S(\d+)_L(\d+)_([RI][A12])_(\d+).fastq(.gz)?$").unwrap();
    static ref BCL2FASTQ_NO_LANE_SPLIT_REGEX: Regex =
        Regex::new(r"^([\w_-]+)_S(\d+)_([RI][A12])_(\d+).fastq(.gz)?$").unwrap();
}

/// Different ways to specify sample names for the `Bcl2FastqDef`
#[derive(Deserialize, Serialize, Clone, PartialEq, Eq, Debug)]
pub enum SampleNameSpec {
    /// All the samples within the fastq directory
    Any,
    /// Only consider a known set of names
    Names(TxHashSet<String>),
}

/// SampleNameSpec for a single sample name
impl From<&str> for SampleNameSpec {
    fn from(sample_name: &str) -> Self {
        let mut names = TxHashSet::default();
        names.insert(sample_name.to_string());
        SampleNameSpec::Names(names)
    }
}

impl SampleNameSpec {
    pub fn contains(&self, sample_name: &str) -> bool {
        match self {
            SampleNameSpec::Any => true,
            SampleNameSpec::Names(ref names) => names.contains(sample_name),
        }
    }
}

/// A pointer to a set of FASTQ files on disk,
/// using the Illumina `bcl2fastq` naming conventions.
/// The `find_fastqs` method will find FASTQ files
/// of the form `heart_1k_v3_S1_L002_R2_001.fastq.gz`
/// with an optional `.gz` suffix.
#[derive(Deserialize, Serialize, Clone, PartialEq, Eq, Debug)]
pub struct Bcl2FastqDef {
    /// The path where to the demulitplexed FASTQ files
    pub fastq_path: String,

    /// Sample name(s) used for this sample
    pub sample_name_spec: SampleNameSpec,

    /// Lanes to include
    pub lane_spec: LaneSpec,
}

impl FindFastqs for Bcl2FastqDef {
    fn find_fastqs(&self) -> Result<Vec<InputFastqs>, Error> {
        let all_fastqs = find_flowcell_fastqs(&self.fastq_path)?;

        let mut res = Vec::new();
        for (info, fastqs) in all_fastqs {
            if self.sample_name_spec.contains(&info.sample)
                && self.lane_spec.contains(info.lane_mode)
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
    pub lane_mode: LaneMode,
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
                    lane_mode: LaneMode::SingleLane(lane),
                    read,
                    chunk,
                    path: path.as_ref().into(),
                });

                return r;
            }

            // Try out the no lane split version next
            if let Some(cap) = BCL2FASTQ_NO_LANE_SPLIT_REGEX.captures(f) {
                let sample = cap.get(1).unwrap().as_str().to_string();
                let s: usize = cap.get(2).unwrap().as_str().parse().unwrap();
                let read = cap.get(3).unwrap().as_str().to_string();
                let chunk: usize = cap.get(4).unwrap().as_str().parse().unwrap();

                let r = Some(IlmnFastqFile {
                    sample,
                    s,
                    lane_mode: LaneMode::NoLaneSplitting,
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

/// Find all the sets of bcl2fastq FASTQ files present in `path` as well as directories directly
/// underneath `path`. Corresponding R1/R2/I1/I2 files are grouped together and reported in an
/// `InputFastqs`, along with a representative `IlmnFastqFile` struct.
pub fn find_flowcell_fastqs(
    path: impl AsRef<Path>,
) -> Result<Vec<(IlmnFastqFile, InputFastqs)>, Error> {
    let mut res = Vec::new();

    let mut files = get_bcl2fastq_files(&path)?;
    // Collect the files which are within the directories underneath `path`.
    // This typically means `path` corresponds to the project folder in the `mkfastq` outs
    for entry in std::fs::read_dir(&path)? {
        let entry = entry?.path();
        if entry.is_dir() {
            files.extend(get_bcl2fastq_files(entry)?);
        }
    }
    files.sort();

    for (_, files) in &files
        .into_iter()
        .group_by(|(info, _)| (info.sample.clone(), info.s, info.lane_mode, info.chunk))
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
    use pretty_assertions::assert_eq;

    #[test]
    fn test_parse() {
        let filename = "heart_1k_v3_S1_L002_R2_001.fastq.gz";
        let r = IlmnFastqFile::new(filename);

        let expected = IlmnFastqFile {
            path: PathBuf::from(filename.to_string()),
            sample: "heart_1k_v3".to_string(),
            s: 1,
            lane_mode: 2.into(),
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
            lane_mode: 2.into(),
            read: "R2".to_string(),
            chunk: 1,
        };

        assert_eq!(r.unwrap(), expected);
    }

    #[test]
    fn test_parse_no_lane_split() {
        let filename = "test/filenames/test_sample_S1_R1_001.fastq.gz";

        let r = IlmnFastqFile::new(filename);

        let expected = IlmnFastqFile {
            path: PathBuf::from(filename),
            sample: "test_sample".to_string(),
            s: 1,
            lane_mode: LaneMode::NoLaneSplitting,
            read: "R1".to_string(),
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
            sample_name_spec: "Infected".into(),
            lane_spec: LaneSpec::Any,
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
            sample_name_spec: "Infected".into(),
            lane_spec: LaneSpec::Lanes(vec![2].into_iter().collect()),
        };

        let fqs = query.find_fastqs()?;
        assert_eq!(fqs.len(), 1);
        assert_eq!(
            fqs[0].r1,
            "test/filenames/bcl2fastq/Infected_S3_L002_R1_001.fastq"
        );
        Ok(())
    }

    #[test]
    fn test_bcl2fastq_no_lane_split() -> Result<(), Error> {
        let path = "test/filenames/bcl2fastq_no_lane_split";

        let query = Bcl2FastqDef {
            fastq_path: path.to_string(),
            sample_name_spec: "test_sample".into(),
            lane_spec: LaneSpec::Any,
        };

        let fqs = query.find_fastqs()?;
        // NOTE: R3 is not used
        let expected = vec![InputFastqs {
            r1: format!("{}/test_sample_S1_R1_001.fastq.gz", path),
            r2: Some(format!("{}/test_sample_S1_R2_001.fastq.gz", path)),
            i1: Some(format!("{}/test_sample_S1_I1_001.fastq.gz", path)),
            i2: None,
            r1_interleaved: false,
        }];
        assert_eq!(fqs, expected);
        Ok(())
    }
}

// The following tests are based on the tests in tenkit: lib/python/tenkit/test/test_fasta.py
#[cfg(test)]
mod tests_from_tenkit {
    use super::*;
    use pretty_assertions::assert_eq;

    #[test]
    fn test_find_input_fastq_files_bcl2fastq_demult() -> Result<(), Error> {
        let path = "test/filenames/bcl2fastq_2";

        let query = Bcl2FastqDef {
            fastq_path: path.to_string(),
            sample_name_spec: "test_sample".into(),
            lane_spec: LaneSpec::Any,
        };

        let fqs = query.find_fastqs()?;
        // NOTE: R3 is not read
        let expected = vec![InputFastqs {
            r1: format!("{}/test_sample_S1_L001_R1_001.fastq.gz", path),
            r2: Some(format!("{}/test_sample_S1_L001_R2_001.fastq.gz", path)),
            i1: Some(format!("{}/test_sample_S1_L001_I1_001.fastq.gz", path)),
            i2: None,
            r1_interleaved: false,
        }];
        assert_eq!(fqs, expected);
        Ok(())
    }

    #[test]
    fn test_find_input_fastq_files_bc2fastq_demult_project_1() -> Result<(), Error> {
        let path = "test/filenames/project_dir";

        let query = Bcl2FastqDef {
            fastq_path: path.to_string(),
            sample_name_spec: "test_sample".into(),
            lane_spec: LaneSpec::Any,
        };

        let fqs = query.find_fastqs()?;
        let expected = vec![
            InputFastqs {
                r1: format!("{}/s1/test_sample_S1_L001_R1_001.fastq.gz", path),
                r2: Some(format!("{}/s1/test_sample_S1_L001_R2_001.fastq.gz", path)),
                i1: Some(format!("{}/s1/test_sample_S1_L001_I1_001.fastq.gz", path)),
                i2: None,
                r1_interleaved: false,
            },
            InputFastqs {
                r1: format!("{}/s2/test_sample_S2_L001_R1_001.fastq.gz", path),
                r2: Some(format!("{}/s2/test_sample_S2_L001_R2_001.fastq.gz", path)),
                i1: Some(format!("{}/s2/test_sample_S2_L001_I1_001.fastq.gz", path)),
                i2: None,
                r1_interleaved: false,
            },
        ];
        assert_eq!(fqs, expected);
        Ok(())
    }

    #[test]
    fn test_find_input_fastq_files_bc2fastq_demult_project_2() -> Result<(), Error> {
        let path = "test/filenames/project_dir";

        let query = Bcl2FastqDef {
            fastq_path: path.to_string(),
            sample_name_spec: "test_sample2".into(),
            lane_spec: LaneSpec::Any,
        };

        let fqs = query.find_fastqs()?;
        let expected = vec![InputFastqs {
            r1: format!("{}/s1/test_sample2_S1_L001_R1_001.fastq.gz", path),
            r2: Some(format!("{}/s1/test_sample2_S1_L001_R2_001.fastq.gz", path)),
            i1: Some(format!("{}/s1/test_sample2_S1_L001_I1_001.fastq.gz", path)),
            i2: None,
            r1_interleaved: false,
        }];
        assert_eq!(fqs, expected);
        Ok(())
    }

    #[test]
    fn test_sample_name_verification() -> Result<(), Error> {
        let path = "test/filenames/tenkit91";
        for &s in ["test_sample", "test_sample_suffix"].iter() {
            let query = Bcl2FastqDef {
                fastq_path: path.to_string(),
                sample_name_spec: s.into(),
                lane_spec: LaneSpec::Any,
            };
            let fqs = query.find_fastqs()?;
            let expected = vec![InputFastqs {
                r1: format!("{}/{}_S1_L001_R1_001.fastq.gz", path, s),
                r2: Some(format!("{}/{}_S1_L001_R2_001.fastq.gz", path, s)),
                i1: Some(format!("{}/{}_S1_L001_I1_001.fastq.gz", path, s)),
                i2: None,
                r1_interleaved: false,
            }];
            assert_eq!(fqs, expected);
        }
        Ok(())
    }

    #[test]
    fn test_sample_name_any() -> Result<(), Error> {
        let path = "test/filenames/tenkit91";

        let query = Bcl2FastqDef {
            fastq_path: path.to_string(),
            sample_name_spec: SampleNameSpec::Any,
            lane_spec: LaneSpec::Any,
        };
        let fqs = query.find_fastqs()?;
        let expected = vec![
            InputFastqs {
                r1: format!("{}/test_sample_S1_L001_R1_001.fastq.gz", path),
                r2: Some(format!("{}/test_sample_S1_L001_R2_001.fastq.gz", path)),
                i1: Some(format!("{}/test_sample_S1_L001_I1_001.fastq.gz", path)),
                i2: None,
                r1_interleaved: false,
            },
            InputFastqs {
                r1: format!("{}/test_sample_suffix_S1_L001_R1_001.fastq.gz", path),
                r2: Some(format!(
                    "{}/test_sample_suffix_S1_L001_R2_001.fastq.gz",
                    path
                )),
                i1: Some(format!(
                    "{}/test_sample_suffix_S1_L001_I1_001.fastq.gz",
                    path
                )),
                i2: None,
                r1_interleaved: false,
            },
        ];
        assert_eq!(fqs, expected);
        Ok(())
    }
}
