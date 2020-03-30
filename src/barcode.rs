// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

//! Tools for dealing with 10x barcodes. Loading barcode whitelists, check barcodes against the whitelist
//! and correcting incorrect barcodes.

use failure::{format_err, Error};
use metric::{SimpleHistogram, TxHashMap, TxHashSet};
use ordered_float::NotNan;
use serde::{Deserialize, Serialize};
use std::hash::Hash;
use std::io::BufRead;
use std::path::Path;

use crate::utils;
use crate::{Barcode, SSeq};

const BC_MAX_QV: u8 = 66; // This is the illumina quality value

/// Return an iterator over the barcode whitelist
fn read_barcode_whitelist_iter(
    filename: impl AsRef<Path>,
) -> Result<impl Iterator<Item = std::io::Result<(SSeq, Option<SSeq>)>>, Error> {
    let barcodes = utils::open_with_gz(filename)?.lines().map(|line| {
        line.and_then(|line| {
            // could be a translation whitelist, take left only
            let mut iter = line.split_whitespace();
            let lhs = SSeq::new(iter.next().unwrap().as_bytes());
            let rhs = iter.next().map(|x| SSeq::new(x.as_bytes()));
            Ok((lhs, rhs))
        })
    });
    Ok(barcodes)
}

/// Read the barcode whitelist and return a vector of SSeq.
pub fn read_barcode_whitelist_vec(filename: impl AsRef<Path>) -> Result<Vec<SSeq>, Error> {
    let barcodes = read_barcode_whitelist_iter(filename)?
        .map(|x| x.map(|x| x.0))
        .collect::<std::io::Result<_>>()?;
    Ok(barcodes)
}

/// Read the barcode whitelist and return a set of SSeq.
pub fn read_barcode_whitelist_set(filename: impl AsRef<Path>) -> Result<TxHashSet<SSeq>, Error> {
    let barcodes = read_barcode_whitelist_iter(filename)?
        .map(|x| x.map(|x| x.0))
        .collect::<std::io::Result<_>>()?;
    Ok(barcodes)
}

/// Read the barcode whitelist and return a map of SSeq to integers.
pub fn read_barcode_whitelist_map(
    filename: impl AsRef<Path>,
) -> Result<TxHashMap<SSeq, u32>, Error> {
    let barcodes = read_barcode_whitelist_iter(filename)?
        .enumerate()
        .map(|(i, x)| x.map(|x| (x.0, i as u32)))
        .collect::<std::io::Result<_>>()?;
    Ok(barcodes)
}

/// Read the barcode (translation) whitelist and return a map of SSeq to SSeq.
pub fn read_barcode_whitelist_trans(
    filename: impl AsRef<Path>,
) -> Result<TxHashMap<SSeq, SSeq>, Error> {
    let barcodes = read_barcode_whitelist_iter(filename.as_ref())?
        .map(|x| {
            let x = x?;
            Ok((
                x.0,
                x.1.ok_or(format_err!(
                    "not a translation whitelist: {:?}",
                    filename.as_ref()
                ))?,
            ))
        })
        .collect::<Result<_, Error>>()?;
    Ok(barcodes)
}

/// Read the barcode whitelist and return a map of SSeq to integers.
#[deprecated = "Use read_barcode_whitelist_map instead"]
pub fn load_barcode_whitelist(filename: impl AsRef<Path>) -> Result<TxHashMap<SSeq, u32>, Error> {
    read_barcode_whitelist_map(filename)
}

pub(crate) fn reduce_counts<K: Hash + Eq>(
    mut v1: TxHashMap<K, u32>,
    mut v2: TxHashMap<K, u32>,
) -> TxHashMap<K, u32> {
    for (k, v) in v2.drain() {
        let counter = v1.entry(k).or_insert(0);
        *counter += v;
    }

    v1
}

#[derive(Serialize, Deserialize)]
pub enum Whitelist {
    Plain(TxHashSet<SSeq>),
    Trans(TxHashMap<SSeq, SSeq>),
}

impl Whitelist {
    pub fn new(p: impl AsRef<Path>) -> Result<Self, Error> {
        // if the parent directory is "translation", we're
        if p.as_ref()
            .parent()
            .and_then(|p| p.file_name().map(|d| d == "translation"))
            .unwrap_or(false)
        {
            Ok(Whitelist::Trans(read_barcode_whitelist_trans(p)?))
        } else {
            Ok(Whitelist::Plain(read_barcode_whitelist_set(p)?))
        }
    }

    fn check(&self, barcode: &mut Barcode, do_translate: bool) -> bool {
        match self {
            Whitelist::Plain(ref whitelist) => {
                if whitelist.contains(&barcode.sequence) {
                    barcode.valid = true;
                } else {
                    barcode.valid = false;
                }
            }
            Whitelist::Trans(ref whitelist) => {
                if let Some(translation) = whitelist.get(&barcode.sequence) {
                    barcode.valid = true;
                    if do_translate {
                        barcode.sequence = *translation;
                    }
                } else {
                    barcode.valid = false;
                }
            }
        }
        barcode.valid
    }
}

/// Check an observed barcode sequence against a whitelist of barcodes
pub struct BarcodeChecker {
    whitelist: Whitelist,
}

impl BarcodeChecker {
    pub fn new(whitelist: impl AsRef<Path>) -> Result<BarcodeChecker, Error> {
        Ok(BarcodeChecker {
            whitelist: Whitelist::new(whitelist)?,
        })
    }

    /// Convert `barcode` to be `is_valid() == true` if the
    /// barcode sequence exactly matches a whitelist sequence.
    pub fn check(&self, barcode: &mut Barcode, do_translate: bool) -> bool {
        self.whitelist.check(barcode, do_translate)
    }
}

/// Load barcode count data from a set of serialized files containing HashMap<Barcode, u32> counts.
pub fn load_barcode_counts(
    count_files: &[impl AsRef<Path>],
) -> Result<TxHashMap<Barcode, u32>, Error> {
    let mut counts = TxHashMap::default();
    for f in count_files {
        let shard_counts = crate::utils::read_obj(f.as_ref())?;
        counts = reduce_counts(counts, shard_counts);
    }

    Ok(counts)
}

const BASE_OPTS: [u8; 4] = [b'A', b'C', b'G', b'T'];
type Of64 = NotNan<f64>;

/// Implement the standard 10x barcode correction algorithm.
/// Requires the barcode whitelist (`whitelist`), and the observed counts
/// of each whitelist barcode, prior to correction (`bc_counts`).
pub struct BarcodeCorrector {
    whitelist: Whitelist,
    bc_counts: SimpleHistogram<Barcode>,
    max_expected_barcode_errors: f64,
    bc_confidence_threshold: f64,
}

impl BarcodeCorrector {
    /// Load a barcode corrector from a whitelist path
    /// and a set of barcode count files.
    /// # Arguments
    /// * `whitelist` - path the barcode whitelist file
    /// * `count_files` - paths to barcode count files written by sort step
    /// * `maximum_expected_barcode_errors` - threshold for sum of probability
    ///   of error on barcode QVs. Reads exceeding this threshold will be marked
    ///   as not valid.
    /// * `bc_confidence_threshold` - if the posterior probability of a correction
    ///    exceeds this threshold, the barcode will be corrected.
    pub fn new(
        whitelist: Whitelist,
        bc_counts: SimpleHistogram<Barcode>,
        max_expected_barcode_errors: f64,
        bc_confidence_threshold: f64,
    ) -> Result<BarcodeCorrector, Error> {
        Ok(BarcodeCorrector {
            whitelist,
            bc_counts,
            max_expected_barcode_errors,
            bc_confidence_threshold,
        })
    }

    /// Attempt to correct a non-whitelist barcode sequence onto the whitelist.
    /// Compute a posterior distribution over whitelist barcodes given:
    /// # A prior distribution over whitelist bacodes (self.counts)
    /// # The likelihood of P(orig_barcode | wl_bc; qual) of the `observed_barcode` sequence
    ///   given a hidden whitelist barcode wl_bc, and the sequencer reported quality values qual.
    /// Use these components to compute a posterior distribution of wl_bc candidates, and correct
    /// the barcode if there's a whitelist barcode with posterior probability greater than
    /// self.bc_confidence_threshold
    #[allow(clippy::needless_range_loop)]
    pub fn correct_barcode(
        &self,
        observed_barcode: &Barcode,
        qual: &[u8],
        do_translate: bool,
    ) -> Option<Barcode> {
        let mut a = observed_barcode.sequence; // Create a copy

        let mut candidates: Vec<(Of64, Barcode)> = Vec::new();
        let mut total_likelihood = Of64::from(0.0);

        for pos in 0..observed_barcode.sequence.len() {
            let qv = qual[pos];
            let existing = a.sequence[pos];
            for val in BASE_OPTS.iter() {
                if *val == existing {
                    continue;
                }
                a.sequence[pos] = *val;
                let mut trial_bc = Barcode {
                    valid: false,
                    sequence: a,
                    gem_group: observed_barcode.gem_group,
                };

                if self.whitelist.check(&mut trial_bc, do_translate) {
                    // Apply additive (Laplace) smoothing.
                    let raw_count = self.bc_counts.get(&trial_bc);
                    let bc_count = 1 + raw_count;
                    let prob_edit = Of64::from(probability(qv.min(BC_MAX_QV)));
                    let likelihood = prob_edit * Of64::from(bc_count as f64);
                    candidates.push((likelihood, trial_bc));
                    total_likelihood += likelihood;
                }
            }
            a.sequence[pos] = existing;
        }

        let thresh = Of64::from(self.bc_confidence_threshold);
        let best_option = candidates.into_iter().max();

        let expected_errors: f64 = qual.iter().cloned().map(probability).sum();

        if let Some((best_like, best_bc)) = best_option {
            if expected_errors < self.max_expected_barcode_errors
                && best_like / total_likelihood >= thresh
            {
                return Some(best_bc);
            }
        }
        None
    }
}

pub fn probability(qual: u8) -> f64 {
    //33 is the illumina qual offset
    let q = f64::from(qual);
    (10_f64).powf(-(q - 33.0) / 10.0)
}

#[cfg(test)]
mod test {
    use super::*;
    use metric::{Metric, SimpleHistogram};
    use proptest::proptest;

    #[test]
    pub fn test_barcode_correction() {
        let mut wl = TxHashSet::default();

        let b1 = Barcode::test_valid(b"AAAAA");
        let b2 = Barcode::test_valid(b"AAGAC");
        let b3 = Barcode::test_valid(b"ACGAA");
        let b4 = Barcode::test_valid(b"ACGTT");

        wl.insert(b1.sequence);
        wl.insert(b2.sequence);
        wl.insert(b3.sequence);
        wl.insert(b4.sequence);

        let mut bc_counts = SimpleHistogram::new();
        bc_counts.insert(b1, 100);
        bc_counts.insert(b2, 11);
        bc_counts.insert(b3, 2);

        let val = BarcodeCorrector {
            max_expected_barcode_errors: 1.0,
            bc_confidence_threshold: 0.95,
            whitelist: Whitelist::Plain(wl),
            bc_counts,
        };

        // Easy
        let t1 = Barcode::test_invalid(b"AAAAA");
        //assert_eq!(val.correct_barcode(&t1, &vec![66,66,66,66,66]), Some(Barcode::test_valid(b"AAAAA")));

        // Low quality
        assert_eq!(
            val.correct_barcode(&t1, &vec![34, 34, 34, 66, 66], false),
            None
        );

        // Trivial correction
        let t2 = Barcode::test_invalid(b"AAAAT");
        assert_eq!(
            val.correct_barcode(&t2, &vec![66, 66, 66, 66, 40], false),
            Some(b1)
        );

        // Pseudo-count kills you
        let t3 = Barcode::test_invalid(b"ACGAT");
        assert_eq!(
            val.correct_barcode(&t3, &vec![66, 66, 66, 66, 66], false),
            None
        );

        // Quality help you
        let t4 = Barcode::test_invalid(b"ACGAT");
        assert_eq!(
            val.correct_barcode(&t4, &vec![66, 66, 66, 66, 40], false),
            Some(b3)
        );

        // Counts help you
        let t5 = Barcode::test_invalid(b"ACAAA");
        assert_eq!(
            val.correct_barcode(&t5, &vec![66, 66, 66, 66, 40], false),
            Some(b1)
        );
    }

    #[test]
    pub fn test_barcode_correction_no_valid_counts() {
        let mut wl = TxHashSet::default();

        let b1 = Barcode::test_valid(b"AAAAA");
        let b2 = Barcode::test_valid(b"AAGAC");
        let b3 = Barcode::test_valid(b"ACGAA");
        let b4 = Barcode::test_valid(b"ACGTT");

        wl.insert(b1.sequence);
        wl.insert(b2.sequence);
        wl.insert(b3.sequence);
        wl.insert(b4.sequence);

        let bc_counts = SimpleHistogram::<Barcode>::new();

        let val = BarcodeCorrector {
            max_expected_barcode_errors: 1.0,
            bc_confidence_threshold: 0.95,
            whitelist: Whitelist::Plain(wl),
            bc_counts,
        };

        // Easy
        let t1 = Barcode::test_invalid(b"AAAAA");
        //assert_eq!(val.correct_barcode(&t1, &vec![66,66,66,66,66]), Some(Barcode::test_valid(b"AAAAA")));

        // Low quality
        assert_eq!(
            val.correct_barcode(&t1, &vec![34, 34, 34, 66, 66], false),
            None
        );

        // Trivial correction
        let t2 = Barcode::test_invalid(b"AAAAT");
        assert_eq!(
            val.correct_barcode(&t2, &vec![66, 66, 66, 66, 40], false),
            Some(b1)
        );
    }

    proptest! {
        #[test]
        fn prop_test_n_in_barcode(
            n_pos in (0..16usize)
        ) {
            let mut wl = TxHashSet::default();
            let bc = Barcode::test_valid(b"GCGATTGACCCAAAGG");
            wl.insert(bc.sequence);

            let corrector = BarcodeCorrector {
                max_expected_barcode_errors: 1.0,
                bc_confidence_threshold: 0.975,
                whitelist: Whitelist::Plain(wl),
                bc_counts: SimpleHistogram::new(),
            };

            let mut bc_seq_with_n = bc.sequence().to_vec();
            bc_seq_with_n[n_pos] = b'N';
            let mut qual = vec![53; bc_seq_with_n.len()];
            qual[n_pos] = 35;
            let bc_with_n = Barcode::test_invalid(&bc_seq_with_n);

            assert_eq!(
                corrector.correct_barcode(&bc_with_n, &qual, false),
                Some(bc)
            );

        }
    }
}
