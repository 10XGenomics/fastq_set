// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

//! Sized, stack-allocated container for a short DNA sequence.

use std::borrow::Borrow;
use std::hash::{Hash, Hasher};
use std::str;

/// Fixed-sized container for a short DNA sequence, up to 23bp in length.
/// Used as a convenient container for barcode or UMI sequences.
#[derive(Serialize, Deserialize, Clone, Copy, PartialOrd, Ord, Eq)]
pub struct SSeq {
    pub(crate) length: u8,
    pub(crate) sequence: [u8; 23],
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
        self.length==0
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
