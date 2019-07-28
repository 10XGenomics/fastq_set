// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

//! Sized, stack-allocated container for a short DNA sequence.

use std::borrow::Borrow;
use std::hash::{Hash, Hasher};
use std::ops::{Index, IndexMut};
use std::str;

use serde::{Serialize, Serializer, Deserialize, Deserializer};
use serde::de::{self, Visitor};

/// Fixed-sized container for a short DNA sequence, up to 23bp in length.
/// Used as a convenient container for barcode or UMI sequences.
#[derive(Clone, Copy, PartialOrd, Ord, Eq)]
pub struct SSeq {
    pub(crate) sequence: [u8; 23],
    pub(crate) length: u8,
}

impl SSeq {
    /// Create a new SSeq from the given byte slice
    pub fn new(seq: &[u8]) -> SSeq {
        assert!(seq.len() <= 23);

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
    pub fn len(&self) -> usize {
        self.length as usize
    }

    pub fn is_empty(&self) -> bool {
        self.length == 0
    }

    pub fn encode_2bit_u32(&self) -> u32 {
        let mut res: u32 = 0;
        assert!(self.length < 16);

        for (bit_pos, str_pos) in (0..self.length).rev().enumerate() {
            let byte: u32 = match self.sequence[str_pos as usize] {
                b'a' => 0,
                b'A' => 0,
                b'c' => 1,
                b'C' => 1,
                b'g' => 2,
                b'G' => 2,
                b't' => 3,
                b'T' => 3,
                _ => panic!("non-ACGT sequence"),
            };

            let v = byte << (bit_pos * 2);

            res = res | v;
        }

        res
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
impl fmt::Debug for SSeq {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut s = String::new();
        for pos in 0..self.len() {
            s.push(self.sequence[pos] as char);
        }
        write!(f, "{}", s)
    }
}

impl Serialize for SSeq {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        serializer.serialize_bytes(self.seq())
    }
}

impl<'de> Deserialize<'de> for SSeq {
    fn deserialize<D>(deserializer: D) -> Result<SSeq, D::Error>
    where
        D: Deserializer<'de>,
    {
       deserializer.deserialize_bytes(SSeqVisitor)
    }
}

struct SSeqVisitor;

impl<'de> Visitor<'de> for SSeqVisitor {
    type Value = SSeq;


    fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        formatter.write_str("an integer between -2^31 and 2^31")
    }


    fn visit_bytes<E>(self, value: &[u8]) -> Result<Self::Value, E>
    where
        E: de::Error,
    {
        Ok(SSeq::new(value))
    }
}

#[cfg(test)]
mod sseq_test {
    use bincode;
    use sseq::SSeq;

    #[test]
    fn sort_test() {
        let s1: &[u8] = b"ASDFQ";
        let s2 = b"ASDFQA";
        let s3 = b"ASDF";
        let s4 = b"ASDF";
        let s5 = b"";
        let s6 = b"A";
        let s7 = b"QERQGHBARTQBEWC";
        let s8 = b"AWERQWERTQER";

        let mut seqs = vec![s1, s2, s3, s4, s5, s6, s7, s8];
        let mut sseqs: Vec<SSeq> = seqs.iter().map(|x| SSeq::new(x)).collect();

        seqs.sort();
        sseqs.sort();

        for i in 0..seqs.len() {
            assert_eq!(seqs[i], sseqs[i].seq());
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
        for i in 0 .. 4 {
            let s = SSeq::new(seq);
            sseqs.push(s);
        }

        let mut buf = Vec::new();
        bincode::serialize_into(&mut buf, &sseqs);
        let roundtrip: Vec<SSeq> = bincode::deserialize_from(&buf[..]).unwrap();
        assert_eq!(sseqs, roundtrip);
    }
}
