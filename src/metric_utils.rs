use bio::pattern_matching;
use metric::PercentMetric;

pub const ILLUMINA_QUAL_OFFSET: u8 = 33;

/// Compute the fraction of 'N' bases in the input sequence
/// as a PercentMetric
#[inline(always)]
pub fn frac_n_bases(seq: &[u8]) -> PercentMetric {
    (
        seq.iter().filter(|&n| *n == b'N').count() as i64,
        seq.len() as i64,
    )
        .into()
}

/// Compute the fraction of Q30 bases in the input sequence
/// quality as a PercentMetric
/// Only bases with quality > 2 is counted in the denominator
pub fn frac_q30_bases(qual: &[u8]) -> PercentMetric {
    let mut num = 0;
    let mut den = 0;
    for q in qual.iter() {
        if *q > 2 + ILLUMINA_QUAL_OFFSET {
            den += 1;
            if *q > 30 + ILLUMINA_QUAL_OFFSET {
                num += 1;
            }
        }
    }
    (num, den).into()
}

pub type Pattern = pattern_matching::bndm::BNDM;
pub struct PatternCheck {
    pattern: Pattern,
}

impl PatternCheck {
    pub fn new(pattern_seq: &[u8]) -> Self {
        PatternCheck {
            pattern: Pattern::new(pattern_seq),
        }
    }
    pub fn exists(&self, read: &[u8]) -> bool {
        self.pattern.find_all(read).next().is_some()
    }
}
