// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

//! Sized, stack-allocated container for a short quality string.

use serde::de::{self, Visitor};
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use std::borrow::Borrow;
use std::hash::{Hash, Hasher};
use std::iter::Iterator;
use std::ops::{Index, IndexMut};
use std::str;

/// Ensure that the input byte slice contains only valid quality characters, and panic otherwise.
pub fn ensure_valid_quality(seq: &[u8]) {
    for (i, &c) in seq.iter().enumerate() {
        let q = c as i16 - 33;
        if q < 0 || q >= 42 {
            panic!(
                "Invalid quality value {} ASCII character {} at position {}",
                q, c, i
            );
        }
    }
}

/// Fixed-sized container for a short quality string, up to 23bp in length.
/// Used as a convenient container for a barcode or UMI quality string.
/// An `SQuality` is guaranteed to contain only valid quality characters.
#[derive(Clone, Copy, PartialOrd, Ord, Eq)]
pub struct SQuality {
    pub(crate) bytes: [u8; 23],
    pub(crate) length: u8,
}

impl SQuality {
    /// Create a new SQuality from the given byte slice
    /// The byte slice must contain only valid quality characters and panics otherwise.
    pub fn new(s: &[u8]) -> SQuality {
        assert!(s.len() <= 23);
        ensure_valid_quality(s);

        let mut bytes = [0u8; 23];
        bytes[0..s.len()].copy_from_slice(&s);

        SQuality {
            length: s.len() as u8,
            bytes,
        }
    }

    /// Returns a byte slice of self.
    pub fn as_bytes(&self) -> &[u8] {
        &self.bytes[0..self.length as usize]
    }

    /// Returns the length of self.
    pub fn len(self) -> usize {
        self.length as usize
    }

    /// Returns true if self has a length of zero bytes.
    pub fn is_empty(self) -> bool {
        self.length == 0
    }
}

impl Index<usize> for SQuality {
    type Output = u8;

    fn index(&self, index: usize) -> &Self::Output {
        if index >= self.length as usize {
            panic!("index out of bounds")
        }

        &self.bytes[index]
    }
}

impl IndexMut<usize> for SQuality {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        if index >= self.length as usize {
            panic!("index out of bounds")
        }

        &mut self.bytes[index]
    }
}

impl AsRef<[u8]> for SQuality {
    fn as_ref(&self) -> &[u8] {
        self.as_bytes()
    }
}

impl Into<String> for SQuality {
    fn into(self) -> String {
        String::from(str::from_utf8(self.as_bytes()).unwrap())
    }
}

impl Borrow<[u8]> for SQuality {
    fn borrow(&self) -> &[u8] {
        self.as_bytes()
    }
}

impl Hash for SQuality {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.as_bytes().hash(state);
    }
}

impl PartialEq for SQuality {
    fn eq(&self, other: &SQuality) -> bool {
        self.as_bytes() == other.as_bytes()
    }
}

use std::fmt;

impl fmt::Display for SQuality {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(str::from_utf8(self.as_bytes()).unwrap())
    }
}
impl fmt::Debug for SQuality {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(&self, f)
    }
}

impl Serialize for SQuality {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        serializer.serialize_bytes(self.as_bytes())
    }
}

impl<'de> Deserialize<'de> for SQuality {
    fn deserialize<D>(deserializer: D) -> Result<SQuality, D::Error>
    where
        D: Deserializer<'de>,
    {
        deserializer.deserialize_bytes(SQualityVisitor)
    }
}

struct SQualityVisitor;

impl<'de> Visitor<'de> for SQualityVisitor {
    type Value = SQuality;

    fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        formatter.write_str("an integer between -2^31 and 2^31")
    }

    fn visit_bytes<E>(self, value: &[u8]) -> Result<Self::Value, E>
    where
        E: de::Error,
    {
        Ok(SQuality::new(value))
    }
}

#[cfg(test)]
mod squality_test {
    use super::*;
    use bincode;
    use proptest::{prop_assert_eq, proptest};

    const VALID_CHARS: &[u8; 42] = br##"!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ"##;

    #[test]
    fn test_sseq_valid_quality() {
        assert_eq!(SQuality::new(&VALID_CHARS[0..23]).len(), 23);
        assert_eq!(SQuality::new(&VALID_CHARS[23..42]).len(), 19);
        assert_eq!(
            SQuality::new(&VALID_CHARS[0..23]).to_string(),
            str::from_utf8(&VALID_CHARS[0..23]).unwrap()
        );
    }

    #[test]
    #[should_panic]
    fn test_sseq_invalid_quality_1() {
        let _ = SQuality::new(b"GHIJ ");
    }

    #[test]
    #[should_panic]
    fn test_sseq_invalid_quality_2() {
        let _ = SQuality::new(b"GHIJK");
    }

    #[test]
    fn test_serde() {
        let mut sseqs = Vec::new();
        sseqs.push(SQuality::new(&VALID_CHARS[0..23]));
        sseqs.push(SQuality::new(&VALID_CHARS[23..41]));

        let mut buf = Vec::new();
        bincode::serialize_into(&mut buf, &sseqs).unwrap();
        let roundtrip: Vec<SQuality> = bincode::deserialize_from(&buf[..]).unwrap();
        assert_eq!(sseqs, roundtrip);
    }

    proptest! {
        #[test]
        fn prop_test_serde_squality(
            ref seq in "[!FGHIJ]{0, 23}",
        ) {
            let target = SQuality::new(seq.as_bytes());
            let encoded: Vec<u8> = bincode::serialize(&target).unwrap();
            let decoded: SQuality = bincode::deserialize(&encoded[..]).unwrap();
            prop_assert_eq!(target, decoded);
        }
    }
}
