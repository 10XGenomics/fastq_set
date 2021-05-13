use serde::de::{self, Visitor};
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use std::borrow::Borrow;
use std::fmt;
use std::hash::{Hash, Hasher};
use std::marker::PhantomData;
use std::ops::{Index, IndexMut};

pub trait ArrayContent {
    fn validate_bytes(bytes: &[u8]);
    fn expected_contents() -> &'static str;
}

/// Fixed-sized container for a short DNA sequence or quality.
/// The capacity is determined by the type `N` and the contents are validates based on type `T`
/// Typically used as a convenient container for barcode or UMI sequences or quality.
#[derive(Clone, Copy, PartialOrd, Ord, Eq)]
pub struct ByteArray<T, const N: usize>
where
    T: ArrayContent,
{
    bytes: [u8; N],
    length: u8,
    phantom: PhantomData<T>,
}

impl<T, const N: usize> ByteArray<T, N>
where
    T: ArrayContent,
{
    pub fn new() -> Self {
        ByteArray {
            length: 0,
            bytes: [0; N],
            phantom: PhantomData,
        }
    }

    /// Caller needs to ensure that the bytes are valid
    pub fn push_unchecked(&mut self, src: &[u8]) {
        let len = self.length as usize;
        assert!(src.len() <= (N - len), "Input slice has length {} which exceeds the remaining capacity of {} bytes in the ByteArray", src.len(), N-len);
        self.bytes[len..len + src.len()].copy_from_slice(&src);
        self.length += src.len() as u8;
    }

    pub fn push(&mut self, src: &[u8]) {
        T::validate_bytes(src);
        self.push_unchecked(src);
    }

    /// Create a new ByteArray from the given byte slice
    /// The byte slice should contain only valid alphabets as defined by ArrayContent trait
    /// otherwise this function will panic
    pub fn from_bytes(src: &[u8]) -> Self {
        let mut arr = Self::new();
        arr.push(src);
        arr
    }

    /// Create a new ByteArray from the given byte slice
    /// Caller needs to ensure that the byte slice contains only valid alphabets as defined by ArrayContent trait
    pub fn from_bytes_unchecked(src: &[u8]) -> Self {
        let mut arr = Self::new();
        arr.push_unchecked(src);
        arr
    }

    pub fn from_iter<'a, C, D>(src: D) -> Self
    where
        C: Borrow<u8>,
        D: IntoIterator<Item = C>,
    {
        let array = ByteArray::from_iter_unchecked(src);
        T::validate_bytes(array.as_bytes());
        array
    }

    pub fn from_iter_unchecked<'a, C, D>(src: D) -> Self
    where
        C: Borrow<u8>,
        D: IntoIterator<Item = C>,
    {
        let mut src = src.into_iter().fuse();
        let mut bytes = [0; N];
        let mut len = 0;
        for (l, r) in bytes.iter_mut().zip(&mut src) {
            *l = *r.borrow();
            len += 1;
        }
        if src.next().is_some() {
            panic!(
                "Error: Input iter exceeds capacity of {} bytes.",
                bytes.len()
            );
        }
        let array = ByteArray {
            length: len,
            bytes,
            phantom: PhantomData,
        };
        array
    }

    /// Returns a byte slice of the contents.
    pub fn as_bytes(&self) -> &[u8] {
        &self.bytes[0..self.length as usize]
    }

    /// Returns a str of the contents.
    pub fn as_str(&self) -> &str {
        std::str::from_utf8(self.as_bytes()).unwrap()
    }

    /// Returns a mutable byte slice of the contents.
    pub fn as_mut_bytes(&mut self) -> &mut [u8] {
        &mut self.bytes[0..self.length as usize]
    }

    /// Returns the length of this sequence, in bytes.
    pub fn len(self) -> usize {
        self.length as usize
    }

    /// Returns true if self has a length of zero bytes.
    pub fn is_empty(self) -> bool {
        self.length == 0
    }

    /// Returns an iterator over the bytes.
    pub fn iter(&self) -> std::slice::Iter<u8> {
        self.as_bytes().iter()
    }
}

impl<T, const N: usize> fmt::Display for ByteArray<T, N>
where
    T: ArrayContent,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_str())
    }
}

impl<T, const N: usize> fmt::Debug for ByteArray<T, N>
where
    T: ArrayContent,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(&self, f)
    }
}

impl<T, const N: usize> Index<usize> for ByteArray<T, N>
where
    T: ArrayContent,
{
    type Output = u8;

    fn index(&self, index: usize) -> &Self::Output {
        if index >= self.length as usize {
            panic!("index out of bounds")
        }

        &self.bytes[index]
    }
}

impl<T, const N: usize> IndexMut<usize> for ByteArray<T, N>
where
    T: ArrayContent,
{
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        if index >= self.length as usize {
            panic!("index out of bounds")
        }
        &mut self.bytes[index]
    }
}

impl<T, const N: usize> AsRef<[u8]> for ByteArray<T, N>
where
    T: ArrayContent,
{
    fn as_ref(&self) -> &[u8] {
        self.as_bytes()
    }
}

impl<T, const N: usize> Into<String> for ByteArray<T, N>
where
    T: ArrayContent,
{
    fn into(self) -> String {
        String::from(self.as_str())
    }
}

impl<T, const N: usize> Borrow<[u8]> for ByteArray<T, N>
where
    T: ArrayContent,
{
    fn borrow(&self) -> &[u8] {
        self.as_bytes()
    }
}

impl<T, const N: usize> Hash for ByteArray<T, N>
where
    T: ArrayContent,
{
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.as_bytes().hash(state);
    }
}

impl<T, const N: usize> PartialEq for ByteArray<T, N>
where
    T: ArrayContent,
{
    fn eq(&self, other: &Self) -> bool {
        self.as_bytes() == other.as_bytes()
    }
}

impl<T, const N: usize> IntoIterator for ByteArray<T, N>
where
    T: ArrayContent,
{
    type Item = u8;
    type IntoIter = std::array::IntoIter<u8, N>;

    fn into_iter(self) -> Self::IntoIter {
        std::array::IntoIter::new(self.bytes)
    }
}

impl<T, const N: usize> Serialize for ByteArray<T, N>
where
    T: ArrayContent,
{
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        serializer.serialize_str(self.as_str())
    }
}

impl<'de, T, const N: usize> Deserialize<'de> for ByteArray<T, N>
where
    T: ArrayContent,
{
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        deserializer.deserialize_str(ByteArrayVisitor {
            phantom_t: PhantomData,
        })
    }
}

struct ByteArrayVisitor<T, const N: usize> {
    phantom_t: PhantomData<[T; N]>,
}

impl<'de, T, const N: usize> Visitor<'de> for ByteArrayVisitor<T, N>
where
    T: ArrayContent,
{
    type Value = ByteArray<T, N>;

    fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        formatter.write_str(T::expected_contents())
    }

    fn visit_str<E>(self, value: &str) -> Result<Self::Value, E>
    where
        E: de::Error,
    {
        Ok(ByteArray::from_bytes(value.as_bytes()))
    }
}
