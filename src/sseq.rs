// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

//! Sized, stack-allocated container for a short DNA sequence.

use serde::de::{self, Visitor};
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use std::borrow::Borrow;
use std::hash::{Hash, Hasher};
use std::iter::Iterator;
use std::ops::{Index, IndexMut};
use std::str;

const UPPER_ACGTN: &[u8; 5] = b"ACGTN";
const N_BASE_INDEX: usize = 4;

/// Make sure that the input byte slice contains only
/// "ACGTN" characters. Panics otherwise with an error
/// message describing the position of the first character
/// that is not an ACGTN.
pub fn ensure_upper_case_acgtn(seq: &[u8]) {
    for (i, &s) in seq.iter().enumerate() {
        if !UPPER_ACGTN.iter().any(|&c| c == s) {
            panic!("Non ACGTN character {} at position {}", s, i);
        };
    }
}

/// Fixed-sized container for a short DNA sequence, up to 23bp in length.
/// Used as a convenient container for barcode or UMI sequences.
/// An `SSeq` is guaranteed to contain only "ACGTN" alphabets
#[derive(Clone, Copy, PartialOrd, Ord, Eq)]
pub struct SSeq {
    pub(crate) sequence: [u8; 23],
    pub(crate) length: u8,
}

impl SSeq {
    /// Create a new SSeq from the given byte slice
    /// The byte slice should contain only "ACGTN" (upper case) alphabets,
    /// otherwise this function will panic
    pub fn new(seq: &[u8]) -> SSeq {
        assert!(seq.len() <= 23);
        ensure_upper_case_acgtn(seq);

        let mut sequence = [0u8; 23];
        sequence[0..seq.len()].copy_from_slice(&seq);

        SSeq {
            length: seq.len() as u8,
            sequence,
        }
    }

    /// Access the sequence data
    pub fn seq(&self) -> &[u8] {
        &self.sequence[0..self.length as usize]
    }

    /// The length of the sequence
    pub fn len(self) -> usize {
        self.length as usize
    }

    pub fn is_empty(self) -> bool {
        self.length == 0
    }

    pub fn encode_2bit_u32(self) -> u32 {
        let mut res: u32 = 0;
        assert!(self.length < 16);

        for (bit_pos, str_pos) in (0..self.length).rev().enumerate() {
            let byte: u32 = match self.sequence[str_pos as usize] {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => panic!("non-ACGT sequence"),
            };

            let v = byte << (bit_pos * 2);

            res |= v;
        }

        res
    }
    pub fn one_hamming_iter(self, opt: HammingIterOpt) -> SSeqOneHammingIter {
        SSeqOneHammingIter::new(self, opt)
    }
}

#[derive(Copy, Clone)]
pub enum HammingIterOpt {
    MutateNBase,
    SkipNBase,
}

/// An iterator over sequences which are one hamming distance away
/// from an `SSeq`. `SSeq` is guaranteed to contain "ACGTN" alphabets.
/// Positions containing "N" or "n" are mutated or skipped
/// depending on the `HammingIterOpt`
pub struct SSeqOneHammingIter {
    source: SSeq,            // Original SSeq from which we need to generate values
    chars: &'static [u8; 5], // Whether it's ACGTN or acgtn
    position: usize,         // Index into SSeq where last base was mutated
    chars_index: usize,      // The last base which was used
    skip_n: bool,            // Whether to skip N bases or mutate them
}

impl SSeqOneHammingIter {
    fn new(sseq: SSeq, opt: HammingIterOpt) -> Self {
        let chars = UPPER_ACGTN;
        SSeqOneHammingIter {
            source: sseq,
            chars,
            position: 0,
            chars_index: 0,
            skip_n: match opt {
                HammingIterOpt::SkipNBase => true,
                HammingIterOpt::MutateNBase => false,
            },
        }
    }
}

impl Iterator for SSeqOneHammingIter {
    type Item = SSeq;

    fn next(&mut self) -> Option<Self::Item> {
        if self.position >= self.source.len() {
            return None;
        }
        let base_at_pos = self.source[self.position];
        if (self.skip_n && base_at_pos == self.chars[N_BASE_INDEX])
            || (self.chars_index >= N_BASE_INDEX)
        {
            // this is an "N" or we went through the ACGT bases already at this position
            self.position += 1;
            self.chars_index = 0;
            self.next()
        } else if base_at_pos == self.chars[self.chars_index] {
            self.chars_index += 1;
            self.next()
        } else {
            let mut next_sseq = self.source;
            next_sseq[self.position] = self.chars[self.chars_index];
            self.chars_index += 1;
            Some(next_sseq)
        }
    }
}

impl Index<usize> for SSeq {
    type Output = u8;

    fn index(&self, index: usize) -> &Self::Output {
        if index >= self.length as usize {
            panic!("index out of bounds")
        }

        &self.sequence[index]
    }
}

impl IndexMut<usize> for SSeq {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        if index >= self.length as usize {
            panic!("index out of bounds")
        }

        &mut self.sequence[index]
    }
}

impl AsRef<[u8]> for SSeq {
    fn as_ref(&self) -> &[u8] {
        self.seq()
    }
}

impl Into<String> for SSeq {
    fn into(self) -> String {
        String::from(str::from_utf8(self.seq()).unwrap())
    }
}

impl Borrow<[u8]> for SSeq {
    fn borrow(&self) -> &[u8] {
        self.seq()
    }
}

impl Hash for SSeq {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.seq().hash(state);
    }
}

impl PartialEq for SSeq {
    fn eq(&self, other: &SSeq) -> bool {
        self.seq() == other.seq()
    }
}

use std::fmt;

impl fmt::Display for SSeq {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut s = String::new();
        for pos in 0..self.len() {
            s.push(self.sequence[pos] as char);
        }
        write!(f, "{}", s)
    }
}
impl fmt::Debug for SSeq {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(&self, f)
    }
}

impl Serialize for SSeq {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        serializer.serialize_str(unsafe { std::str::from_utf8_unchecked(self.seq()) })
    }
}

impl<'de> Deserialize<'de> for SSeq {
    fn deserialize<D>(deserializer: D) -> Result<SSeq, D::Error>
    where
        D: Deserializer<'de>,
    {
        deserializer.deserialize_str(SSeqVisitor)
    }
}

struct SSeqVisitor;

impl<'de> Visitor<'de> for SSeqVisitor {
    type Value = SSeq;

    fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        formatter.write_str("An [ACGTN]* string")
    }

    fn visit_str<E>(self, value: &str) -> Result<Self::Value, E>
    where
        E: de::Error,
    {
        Ok(SSeq::new(value.as_bytes()))
    }
}

#[cfg(test)]
mod sseq_test {
    use super::*;
    use bincode;
    use itertools::{assert_equal, Itertools};
    use proptest::collection::vec;
    use proptest::{prop_assert_eq, proptest};
    use std::collections::HashSet;

    #[test]
    fn sort_test() {
        let s1: &[u8] = b"ACNGTA";
        let s2 = b"TAGTCGGC";
        let s3 = b"CATC";
        let s4 = b"TGTG";
        let s5 = b"";
        let s6 = b"A";
        let s7 = b"AACCATAGCCGGNATC";
        let s8 = b"GAACNAGNTGGA";

        let mut seqs = vec![s1, s2, s3, s4, s5, s6, s7, s8];
        let mut sseqs: Vec<SSeq> = seqs.iter().map(|x| SSeq::new(x)).collect();

        seqs.sort();
        sseqs.sort();

        for i in 0..seqs.len() {
            assert_eq!(seqs[i], sseqs[i].seq());
        }
    }

    proptest! {
        #[test]
        fn prop_test_sort_sseq(
            ref seqs_str in vec("[ACGTN]{0, 23}", 0usize..=10usize),
        ) {
            let mut seqs = seqs_str.iter().map(|s| s.clone().into_bytes()).collect_vec();
            let mut sseqs: Vec<SSeq> = seqs.iter().map(|x| SSeq::new(x)).collect();

            seqs.sort();
            sseqs.sort();

            for i in 0..seqs.len() {
                assert_eq!(seqs[i], sseqs[i].seq());
            }
        }
    }

    #[test]
    fn dna_encode() {
        let s1 = SSeq::new(b"AAAAA");
        assert_eq!(s1.encode_2bit_u32(), 0);

        let s1 = SSeq::new(b"AAAAT");
        assert_eq!(s1.encode_2bit_u32(), 3);

        let s1 = SSeq::new(b"AAACA");
        assert_eq!(s1.encode_2bit_u32(), 4);

        let s1 = SSeq::new(b"AACAA");
        assert_eq!(s1.encode_2bit_u32(), 16);

        let s1 = SSeq::new(b"AATA");
        assert_eq!(s1.encode_2bit_u32(), 12);
    }

    #[test]
    fn test_serde() {
        let seq = b"AGCTAGTCAGTCAGTA";
        let mut sseqs = Vec::new();
        for _ in 0..4 {
            let s = SSeq::new(seq);
            sseqs.push(s);
        }

        let mut buf = Vec::new();
        bincode::serialize_into(&mut buf, &sseqs).unwrap();
        let roundtrip: Vec<SSeq> = bincode::deserialize_from(&buf[..]).unwrap();
        assert_eq!(sseqs, roundtrip);
    }

    #[test]
    fn test_serde_json() {
        let seq = SSeq::new(b"AGCTAGTCAGTCAGTA");
        let json_str = serde_json::to_string(&seq).unwrap();
        assert_eq!(json_str, r#""AGCTAGTCAGTCAGTA""#);
    }

    proptest! {
        #[test]
        fn prop_test_serde_sseq(
            ref seq in "[ACGTN]{0, 23}",
        ) {
            let target = SSeq::new(seq.as_bytes());
            let encoded: Vec<u8> = bincode::serialize(&target).unwrap();
            let decoded: SSeq = bincode::deserialize(&encoded[..]).unwrap();
            prop_assert_eq!(target, decoded);
        }
        #[test]
        fn prop_test_serde_json_sseq(
            ref seq in "[ACGTN]{0, 23}",
        ) {
            let target = SSeq::new(seq.as_bytes());
            let encoded = serde_json::to_string_pretty(&target).unwrap();
            let decoded: SSeq = serde_json::from_str(&encoded).unwrap();
            prop_assert_eq!(target, decoded);
        }
    }

    fn test_hamming_helper(seq: &String, opt: HammingIterOpt, n: u8) {
        let sseq = SSeq::new(seq.as_bytes());
        // Make sure that the hamming distance is 1 for all elements
        for neighbor in sseq.one_hamming_iter(opt) {
            assert_eq!(
                sseq.seq()
                    .iter()
                    .zip_eq(neighbor.seq().iter())
                    .filter(|(a, b)| a != b)
                    .count(),
                1
            );
        }
        // Make sure that the total number of elements is equal to what we expect.
        let non_n = sseq.seq().iter().filter(|&&s| s != n).count();
        let n_bases = sseq.len() - non_n;
        assert_eq!(
            sseq.one_hamming_iter(opt).collect::<HashSet<_>>().len(),
            match opt {
                HammingIterOpt::SkipNBase => 3 * non_n,
                HammingIterOpt::MutateNBase => 3 * non_n + 4 * n_bases,
            }
        );
    }

    proptest! {
        #[test]
        fn prop_test_one_hamming_upper(
            seq in "[ACGTN]{0, 23}", // 0 and 23 are inclusive bounds
        ) {
            test_hamming_helper(&seq, HammingIterOpt::SkipNBase, b'N');
            test_hamming_helper(&seq, HammingIterOpt::MutateNBase, b'N');
        }
    }

    #[test]
    #[should_panic]
    fn test_sseq_invalid_1() {
        let _ = SSeq::new(b"ASDF");
    }

    #[test]
    #[should_panic]
    fn test_sseq_invalid_2() {
        let _ = SSeq::new(b"ag");
    }

    #[test]
    fn test_one_hamming_simple() {
        assert_equal(
            SSeq::new(b"GAT")
                .one_hamming_iter(HammingIterOpt::SkipNBase)
                .collect_vec(),
            vec![
                b"AAT", b"CAT", b"TAT", b"GCT", b"GGT", b"GTT", b"GAA", b"GAC", b"GAG",
            ]
            .into_iter()
            .map(|x| SSeq::new(x)),
        );

        assert_equal(
            SSeq::new(b"GNT")
                .one_hamming_iter(HammingIterOpt::SkipNBase)
                .collect_vec(),
            vec![b"ANT", b"CNT", b"TNT", b"GNA", b"GNC", b"GNG"]
                .into_iter()
                .map(|x| SSeq::new(x)),
        );
    }
}
