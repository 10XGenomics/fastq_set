/// Fixed-sized container for a short DNA sequence, up to 23bp in length.
/// Used as a convenient container for barcode or UMI sequences.
#[derive(Serialize, Deserialize, Clone, Copy, PartialOrd, Ord, PartialEq, Eq, Hash)]
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
}

impl AsRef<[u8]> for SSeq {
    fn as_ref(&self) -> &[u8] {
        self.seq()
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
