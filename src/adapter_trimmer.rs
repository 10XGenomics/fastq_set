
use bio::alignment::sparse;
use bio::alignment::pairwise::{self, MatchParams, Scoring};
use bio::alignment::Alignment;
use std::cmp::{min, max};
use std::i32;
use bio::alignment::sparse::HashMapFx;
use std::ops::Range;
use WhichEnd;

type Aligner = pairwise::banded::Aligner<MatchParams>;

const KMER_LEN: usize = 5;
const WINDOW_LEN: usize = 2;
const EXPECTED_READ_LEN: usize = 150;
const ALLOWED_ERROR_RATE: f64 = 0.1_f64;
const MATCH_SCORE: i32 = 1;
const EDIT_SCORE: i32 = -2;
const MIN_PATH_SCORE: i32 = -5;

#[derive(Serialize, Deserialize, Debug, Copy, Clone, PartialEq)]
pub enum AdapterLoc {
    #[serde(rename = "anywhere")]
    Anywhere,
    #[serde(rename = "non_internal")]
    NonInternal,
    #[serde(rename = "anchored")]
    Anchored,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct Adapter {
    pub name: String,
    pub end: WhichEnd,
    pub location: AdapterLoc,
    pub seq: String,
}

impl Adapter {
    pub fn new(name: String, end: WhichEnd, location: AdapterLoc, seq: String) -> Self {
        Adapter {
            name,
            end,
            location,
            seq,
        }
    }

    pub fn seq(&self) -> &String {
        &self.seq
    }
}

pub struct CutAdapt<'a> {
    adapter: &'a Adapter,
    seq: &'a [u8],
    kmer_hash: HashMapFx<&'a[u8], Vec<u32>>,
    aligner: Aligner,
    cut_scores: CutScores,
}

fn count_errors(alignment: &Alignment) -> usize {
    use bio::alignment::AlignmentOperation::{Del, Ins, Subst};
    let clip_errors = min(alignment.ylen - alignment.yend, alignment.xlen - alignment.xend);
    clip_errors + alignment.operations.iter().filter(|&op| *op==Subst || *op==Ins || *op==Del).count()
}
fn count_matches(alignment: &Alignment) -> usize {
    use bio::alignment::AlignmentOperation::Match;
    alignment.operations.iter().filter(|&op| *op==Match).count()
}
fn count_ops(alignment: &Alignment) -> usize {
    use bio::alignment::AlignmentOperation::{Del, Ins, Match, Subst};
    alignment.operations.iter().filter(|&op| *op==Match || *op==Subst || *op==Ins || *op==Del).count()
}

impl<'a> CutAdapt<'a> {
    pub fn new(adapter: &'a Adapter) -> Self {

        let cut_scores = CutScores::new(adapter.end, adapter.location);
        let scoring = Scoring::from_scores(0, EDIT_SCORE, MATCH_SCORE, EDIT_SCORE)
            .xclip_prefix(cut_scores.read_prefix)
            .xclip_suffix(cut_scores.read_suffix)
            .yclip_prefix(cut_scores.adapter_prefix)
            .yclip_suffix(cut_scores.adapter_suffix);

        let aligner = Aligner::with_capacity_and_scoring(
            EXPECTED_READ_LEN,
            adapter.seq.len(),
            scoring,
            KMER_LEN,
            WINDOW_LEN,
        );
        CutAdapt {
            adapter,
            seq: adapter.seq.as_bytes(),
            kmer_hash: sparse::hash_kmers(adapter.seq.as_bytes(), KMER_LEN),
            aligner,
            cut_scores,
        }

    }

    pub fn cutadapt(&mut self, read: &[u8]) -> Option<CutAdaptResult> {
        let matches = sparse::find_kmer_matches_seq2_hashed(read, &self.kmer_hash, KMER_LEN);
        let matches = sparse::expand_kmer_matches(read, self.seq, KMER_LEN, &matches, 1);
        let band_path = compute_path(&matches, read.len(), self.seq.len(), KMER_LEN, &self.cut_scores);
        match band_path {
            Some(p) => {
                let score = p.total_score();
                if score > MIN_PATH_SCORE {

                    // Set the clipping score correctly if the
                    // adapter location is `Anywhere`
                    if self.adapter.location == AdapterLoc::Anywhere {
                        let scoring = self.aligner.get_mut_scoring();
                        match self.adapter.end {
                            WhichEnd::ThreePrime => {
                                let dx = read.len() - p.last_match.0 as usize;
                                let dy = self.adapter.seq.len() - p.last_match.1 as usize;
                                // if (max(dx, dy) - min(dx, dy)) < 3 {
                                if dx==dy {
                                    scoring.xclip_suffix = pairwise::MIN_SCORE;
                                    scoring.yclip_suffix = pairwise::MIN_SCORE;
                                } else if dx > dy {
                                    scoring.xclip_suffix = 0;
                                    scoring.yclip_suffix = pairwise::MIN_SCORE;
                                } else {
                                    scoring.xclip_suffix = pairwise::MIN_SCORE;
                                    scoring.yclip_suffix = 0;
                                }
                            },
                            WhichEnd::FivePrime => {
                                let dx = p.first_match.0;
                                let dy = p.first_match.1;
                                if dx==dy {
                                    scoring.xclip_prefix = pairwise::MIN_SCORE;
                                    scoring.yclip_prefix = pairwise::MIN_SCORE;
                                } else if dx > dy {
                                    scoring.xclip_prefix = 0;
                                    scoring.yclip_prefix = pairwise::MIN_SCORE;
                                } else {
                                    scoring.xclip_prefix = pairwise::MIN_SCORE;
                                    scoring.yclip_prefix = 0;
                                }
                            },
                        }
                    }
                    let alignment = self.aligner.custom_with_match_path(read, self.seq, &matches, &p.path);
                    let max_score = min(alignment.ylen, count_ops(&alignment));
                    // let max_score = min(alignment.ylen, alignment.y_aln_len() + (alignment.xlen - alignment.xend));
                    let err_rate = (count_errors(&alignment) as f64) / (max_score as f64);
                    // println!("{} {}", max_score, err_rate);
                    // println!("{}", alignment.pretty(read, &self.seq));
                    // let match_path = p.path.iter().map(|x| matches[*x]).collect::<Vec<_>>();
                    // println!("{:?}", match_path);
                    // self.aligner.visualize(&alignment);
                    // println!("{:?}", alignment);
                    if err_rate <= ALLOWED_ERROR_RATE && count_matches(&alignment) >= KMER_LEN {
                        Some(CutAdaptResult {
                            start: alignment.xstart,
                            end: alignment.xend,
                            path_score: score,
                            trim_range: match self.adapter.end {
                                WhichEnd::ThreePrime => (0..alignment.xstart),
                                WhichEnd::FivePrime => (alignment.xend..alignment.xlen),
                            }
                        })
                    } else {
                        None
                    }
                } else {
                    None
                }
            },
            None => None
        }
    }
}

pub struct CutAdaptResult {
    pub start: usize,
    pub end: usize,
    pub path_score: i32,
    pub trim_range: Range<usize>,
}

#[derive(Debug, Copy, Clone)]
struct CutScores {
    match_score: i32,
    edit_score: i32,
    read_prefix: i32,
    read_suffix: i32,
    adapter_prefix: i32,
    adapter_suffix: i32,
    right_cut: i32,
    left_cut: i32,
}

impl CutScores {
    fn new(end: WhichEnd, loc: AdapterLoc) -> Self {
        use self::WhichEnd::{FivePrime, ThreePrime};
        use self::AdapterLoc::{Anchored, Anywhere, NonInternal};
        let mut scores = CutScores {
            match_score: MATCH_SCORE,
            edit_score: EDIT_SCORE,
            read_prefix: pairwise::MIN_SCORE,
            read_suffix: pairwise::MIN_SCORE,
            adapter_prefix: pairwise::MIN_SCORE,
            adapter_suffix: pairwise::MIN_SCORE,
            right_cut: pairwise::MIN_SCORE,
            left_cut: pairwise::MIN_SCORE,
        };

        match (end, loc) {
            (ThreePrime, Anywhere) => {
                scores.read_prefix = 0;
                scores.right_cut = 0;
            },
            (ThreePrime, NonInternal) => {
                scores.read_prefix = 0;
                scores.adapter_suffix = 0;
            },
            (ThreePrime, Anchored) => {
                scores.read_prefix = 0;
            },
            (FivePrime, Anywhere) => {
                scores.read_suffix = 0;
                scores.left_cut = 0;
            },
            (FivePrime, NonInternal) => {
                scores.read_suffix = 0;
                scores.adapter_prefix = 0;
            },
            (FivePrime, Anchored) => {
                scores.read_suffix = 0;
            }
        }

        scores
    }
}


#[derive(Debug, Clone)]
struct PathTracker {
    parent: Option<usize>,
    index: usize,
    last_match: (u32, u32),
    left_score: i32,
    right_score: i32,
    internal_score: i32,
}
impl PathTracker {
    fn total_score(&self) -> i32 {
        self.left_score + self.internal_score + self.right_score
    }
}

#[derive(Debug, Clone)]
pub struct PathData {
    pub path: Vec<usize>,
    pub first_match: (u32, u32),
    pub last_match: (u32, u32),
    left_score: i32,
    right_score: i32,
    internal_score: i32,
}
impl PathData {
    pub fn total_score(&self) -> i32 {
        self.left_score + self.internal_score + self.right_score
    }
}

fn compute_path(
    matches: &[(u32, u32)],
    read_len: usize,
    adapter_len: usize,
    k: usize,
    cut_scores: &CutScores,
) -> Option<PathData> {

    if matches.is_empty() {
        return None;
    }

    let max_allowed_indel = (adapter_len as f64 * ALLOWED_ERROR_RATE) as u32;

    // incoming matches must be sorted.
    for i in 1..matches.len() {
        assert!(matches[i - 1] < matches[i]);
    }
    // println!("Matches -> {:?}", matches);
    let mut best_path_ending_here: Vec<PathTracker> = Vec::with_capacity(matches.len());

    for i in 0..matches.len() {
        let this_match = matches[i];

        let left_score = {
            let read_score = max((this_match.0 as i32) * cut_scores.edit_score, cut_scores.read_prefix);
            let adapter_score = max((this_match.1 as i32) * cut_scores.edit_score, cut_scores.adapter_prefix);
            max(read_score + adapter_score, cut_scores.left_cut + max(read_score, adapter_score))
        };

        let right_score = {
            let read_edit_score = ((read_len - k - (this_match.0 as usize)) as i32) * cut_scores.edit_score;
            let adapter_edit_score = ((adapter_len - k - (this_match.1 as usize)) as i32) * cut_scores.edit_score;
            let read_score = max(read_edit_score, cut_scores.read_suffix);
            let adapter_score = max(adapter_edit_score, cut_scores.adapter_suffix);
            max(read_score + adapter_score, cut_scores.right_cut + max(read_score, adapter_score))
        };

        let mut best_path = PathTracker {
            parent: None,
            index: i,
            last_match: this_match.clone(),
            left_score,
            right_score,
            internal_score: (k as i32) * cut_scores.match_score,
        };

        // println!("{:?}", this_match);
        for j in 0..i {
            let p = &best_path_ending_here[j];
            if (p.last_match.0+1 == this_match.0) && (p.last_match.1+1 == this_match.1) {
                // println!("Continues");
                let score = p.left_score + p.internal_score + cut_scores.match_score + right_score;
                if score > best_path.total_score() {
                    best_path.parent = Some(j);
                    best_path.left_score = p.left_score;
                    best_path.internal_score += cut_scores.match_score;
                }
            } else if ((p.last_match.0 + k as u32) <= this_match.0) && ((p.last_match.1 + k as u32) <= this_match.1) {
                let dx = this_match.0 - p.last_match.0 - k as u32;
                let dy = this_match.1 - p.last_match.1 - k as u32;
                let gap = max(dx, dy);
                let min_indel = gap - min(dx, dy);
                let internal_score = p.internal_score + (gap as i32) * cut_scores.edit_score + (k as i32) * cut_scores.match_score;
                let score = p.left_score + internal_score + right_score;
                if score > best_path.total_score() && min_indel <= max_allowed_indel {
                    best_path.parent = Some(j);
                    best_path.left_score = p.left_score;
                    best_path.internal_score = internal_score;
                }
            } 
        }
        // println!("{:?}", best_path);
        best_path_ending_here.push(best_path);
    }

    // println!("{:?}", best_path_ending_here);
    let mut best_end = best_path_ending_here.iter().max_by_key(|x| x.total_score()).unwrap();
    let best_ind = best_end.index;

    let mut path = Vec::with_capacity(matches.len());
    path.push(best_end.index);
    while let Some(p) = best_end.parent {
        path.push(p);
        best_end = &best_path_ending_here[p];
    }
    path.reverse();
    
    let first_match = matches[path[0]].clone();
    Some(PathData {
        path,
        first_match,
        last_match: best_path_ending_here[best_ind].last_match,
        left_score: best_path_ending_here[best_ind].left_score,
        right_score: best_path_ending_here[best_ind].right_score,
        internal_score: best_path_ending_here[best_ind].internal_score,
    })
    
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;
    #[test]
    fn test_anywhere_3p() {
        test_helper(WhichEnd::ThreePrime, AdapterLoc::Anywhere, "tests/cutadapt_v1_18_anywhere_3p.fa");
    }

    #[test]
    fn test_anywhere_5p() {
        test_helper(WhichEnd::FivePrime, AdapterLoc::Anywhere, "tests/cutadapt_v1_18_anywhere_5p.fa");
    }

    #[test]
    fn test_anchored_3p() {
        test_helper(WhichEnd::ThreePrime, AdapterLoc::Anchored, "tests/cutadapt_v1_18_anchored_3p.fa");
    }

    #[test]
    fn test_anchored_5p() {
        test_helper(WhichEnd::FivePrime, AdapterLoc::Anchored, "tests/cutadapt_v1_18_anchored_5p.fa");
    }

    #[test]
    fn test_non_internal_3p() {
        test_helper(WhichEnd::ThreePrime, AdapterLoc::NonInternal, "tests/cutadapt_v1_18_non_internal_3p.fa");
    }

    #[test]
    fn test_non_internal_5p() {
        test_helper(WhichEnd::FivePrime, AdapterLoc::NonInternal, "tests/cutadapt_v1_18_non_internal_5p.fa");
    }

    fn test_helper(end: WhichEnd, loc: AdapterLoc, ref_path: impl AsRef<Path>) {
        use bio::io::fasta::Reader;

        let adapter = Adapter::new("primer".into(), end, loc, "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC".into());
        let mut cutadapt = CutAdapt::new(&adapter);

        let input = Reader::from_file("tests/input.fa").unwrap();
        let output = Reader::from_file(ref_path).unwrap();

        let mut total = 0;
        let mut correct_untrimmed = 0;
        let mut self_trimmed_only = 0;
        let mut self_untrimmed_only = 0;
        let mut correct_trimmed = 0;
        let mut inconsistent_trimmed = 0;

        for (r_i, r_o) in input.records().zip(output.records()) {
            let r_i = r_i.unwrap();
            let r_o = r_o.unwrap();

            total += 1;

            let seq = r_i.seq();
            let seq_len = seq.len();
            let expected_trimmed_len = r_o.seq().len();

            let cutadapt_result = cutadapt.cutadapt(seq);
            let trimmed_len = match cutadapt_result {
                Some(r) => r.trim_range.len(),
                None => seq.len(),
            };

            if seq_len==expected_trimmed_len && trimmed_len==seq_len {
                correct_untrimmed += 1;
            } else if seq_len==expected_trimmed_len {
                self_trimmed_only += 1;
            } else if trimmed_len==seq_len {
                self_untrimmed_only += 1;
            } else if (max(expected_trimmed_len, trimmed_len) - min(expected_trimmed_len, trimmed_len)) <= 2 {
                correct_trimmed += 1;
            } else {
                inconsistent_trimmed += 1;
            }

        }

        // Sensitivity = (Self & Ref) / Ref
        let sensitivity = (correct_trimmed as f64 + inconsistent_trimmed as f64) / (correct_trimmed as f64 + inconsistent_trimmed as f64 + self_untrimmed_only as f64);
        // PPV = (Self & Ref)/Self
        let ppv = (correct_trimmed as f64 + inconsistent_trimmed as f64) / (correct_trimmed as f64 + inconsistent_trimmed as f64 + self_trimmed_only as f64);
        // Concordance = (Self & Ref & Self==Ref) / (Self & Ref)
        let concordance = (correct_trimmed as f64) / (correct_trimmed as f64 + inconsistent_trimmed as f64);

        println!("{:?}", adapter);
        println!("Reads not trimmed {}/{}", correct_untrimmed, total);
        println!("Sensitivity {:.2}%", 100.0*sensitivity);
        println!("PPV         {:.2}%", 100.0*ppv);
        println!("Concordance {:.2}%", 100.0*concordance);

        assert!(sensitivity > 0.99_f64);
        assert!(ppv > 0.99_f64);
        assert!(concordance > 0.99_f64);
    }
}