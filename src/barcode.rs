// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

use std::cmp::max;
use ordered_float::NotNaN;
use failure::Error;
use std::path::Path;

use std::io::BufRead;
use fxhash::{FxHashSet, FxHashMap};

use {Barcode, SSeq};


/// Load a (possibly gzipped) barcode whitelist file.
/// Each line in the file is a single whitelist barcode.
fn load_barcode_whitelist(filename: impl AsRef<Path>) -> Result<FxHashMap<SSeq, u32>, Error> {

    let reader = ::fastq_read::open_w_gz(filename)?;
    let mut bc_map = FxHashMap::default();

    let mut i = 1u32;
    for l in reader.lines() {
        let seq = SSeq::new(&l?.into_bytes());
        bc_map.insert(seq, i);
        i += 1;
    }

    Ok(bc_map)
}

pub struct BarcodeChecker {
    whitelist: FxHashSet<SSeq>
}

impl BarcodeChecker { 
    pub fn new(whitelist: impl AsRef<Path>) -> Result<BarcodeChecker, Error> {
        let wl = load_barcode_whitelist(whitelist)?;

        let mut whitelist = FxHashSet::default();
        whitelist.extend(wl.keys());
        Ok(BarcodeChecker { whitelist })
    }

    pub fn check(&self, barcode: &mut Barcode) {
        if self.whitelist.contains(&barcode.sequence) {
            barcode.valid = true;
        } else {
            barcode.valid = false;
        }
    }
}



pub struct BarcodeCorrector<'a> {
	pub whitelist: &'a FxHashMap<Vec<u8>, u32>,
    pub bc_counts: &'a FxHashMap<(u8, u32), u32>,
	pub max_expected_barcode_errors: f64,
	pub bc_confidence_threshold: f64,
}

const BASE_OPTS: [u8; 4] = [b'A', b'C', b'G', b'T'];

type Of64 = NotNaN<f64>;

impl<'a> BarcodeCorrector<'a> {

    fn do_correct_barcode(&self, gem_group: u8, barcode: &[u8], qual: &[u8]) -> Option<u32> {
        let mut a = Vec::from(barcode);

        let mut candidates: Vec<(Of64, Vec<u8>)> = Vec::new();
        let mut total_likelihood = Of64::from(0.0);
         
	    for pos in 0 .. barcode.len() {
            let qv = qual[pos];
            let existing = a[pos];
            for val in BASE_OPTS.iter().cloned() {

                if val == existing { continue; }
                a[pos] = val;

                match self.whitelist.get(&a) {
                    Some(bc_id) => {
                        let bc_count = self.bc_counts.get(&(gem_group, *bc_id)).cloned().unwrap_or(0);
                        let prob_edit = max(Of64::from(0.0005), Of64::from(probability(qv)));
                        let likelihood = prob_edit * max(Of64::from(bc_count as f64), Of64::from(0.5));
                        candidates.push((likelihood, a.clone()));
                        total_likelihood += likelihood;
                    },
                    None => (),
                }
            }
            a[pos] = existing;
        }
	
        println!("{:?}", candidates);

        let thresh = Of64::from(self.bc_confidence_threshold);

        let best_option = candidates.into_iter().max();
        
        match best_option {
            Some((best_like, bc)) => {
                if best_like / total_likelihood > thresh {
                    self.whitelist.get(&bc).cloned()
                } else {
                    None
                }
            },
            _ => None
        }
    }

    pub fn correct_barcode(&self, gem_group: u8, barcode: &[u8], qual: &[u8]) -> Option<u32> {
        let expected_errors: f64 = qual.iter().cloned().map(probability).sum();

        let bc_id = 
            match self.whitelist.get(barcode) {
                Some(id) => Some(*id),
                None => self.do_correct_barcode(gem_group, barcode, qual),
            };

        if bc_id.is_some() && expected_errors < self.max_expected_barcode_errors {
            bc_id
        } else {
            None
        }
    }
}


pub fn probability(qual: u8) -> f64 {
    //33 is the illumina qual offset
    let q = qual as f64;
    (10_f64).powf(-(q - 33.0) / 10.0) 
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    pub fn test_bc_correct()
    {
        let mut wl = FxHashMap::default();
        wl.insert(Vec::from(b"AAAAA".as_ref()), 1);
        wl.insert(Vec::from(b"AAGAC".as_ref()), 2);
        wl.insert(Vec::from(b"ACGAA".as_ref()), 3);
        wl.insert(Vec::from(b"ACGTT".as_ref()), 4);

        let mut counts = FxHashMap::default();
        counts.insert((0,1), 100);
        counts.insert((0,2), 11);
        counts.insert((0,3), 2);

        let val = BarcodeCorrector {
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