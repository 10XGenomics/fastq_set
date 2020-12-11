pub use generic_array::typenum;
use generic_array::{ArrayLength, GenericArray};
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
/// Typically used as a convenient container for barcode or UMI sequences.
#[derive(Clone, Copy, PartialOrd, Ord, Eq)]
pub struct ByteArray<N, T>
where
    N: ArrayLength<u8>,
    N::ArrayType: Copy,
    T: ArrayContent,
{
    bytes: GenericArray<u8, N>,
    length: u8,
    phantom: PhantomData<T>,
}

impl<N, T> ByteArray<N, T>
where
    N: ArrayLength<u8>,
    N::ArrayType: Copy,
    T: ArrayContent,
{
    /// Create a new ByteArray from the given byte slice
    /// The byte slice should contain only valid alphabets as defined by ArrayContent trait
    /// otherwise this function will panic
    pub fn new(src: &[u8]) -> Self {
        let mut bytes: GenericArray<u8, N> = GenericArray::default();

        assert!(src.len() <= bytes.len());
        T::validate_bytes(src);

        for (l, r) in bytes.as_mut_slice().iter_mut().zip(src.iter()) {
            *l = *r;
        }
        ByteArray {
            length: src.len() as u8,
            bytes,
            phantom: PhantomData,
        }
    }

    /// Returns a byte slice of the contents.
    pub fn as_bytes(&self) -> &[u8] {
        &self.bytes.as_slice()[0..self.length as usize]
    }

    /// Returns a str of the contents.
    pub fn as_str(&self) -> &str {
        std::str::from_utf8(self.as_bytes()).unwrap()
    }

    /// Returns a mutable byte slice of the contents.
    pub fn as_mut_bytes(&mut self) -> &mut [u8] {
        &mut self.bytes.as_mut_slice()[0..self.length as usize]
    }

    /// Returns the length of this sequence, in bytes.
    pub fn len(self) -> usize {
        self.length as usize
    }

    /// Returns an iterator over the bytes.
    pub fn iter(&self) -> std::slice::Iter<u8> {
        self.as_bytes().iter()
    }
}

impl<N, T> fmt::Display for ByteArray<N, T>
where
    N: ArrayLength<u8>,
    N::ArrayType: Copy,
    T: ArrayContent,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_str())
    }
}

impl<N, T> fmt::Debug for ByteArray<N, T>
where
    N: ArrayLength<u8>,
    N::ArrayType: Copy,
    T: ArrayContent,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(&self, f)
    }
}

impl<N, T> Index<usize> for ByteArray<N, T>
where
    N: ArrayLength<u8>,
    N::ArrayType: Copy,
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

impl<N, T> IndexMut<usize> for ByteArray<N, T>
where
    N: ArrayLength<u8>,
    N::ArrayType: Copy,
    T: ArrayContent,
{
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        if index >= self.length as usize {
            panic!("index out of bounds")
        }
        &mut self.bytes[index]
    }
}

impl<N, T> AsRef<[u8]> for ByteArray<N, T>
where
    N: ArrayLength<u8>,
    N::ArrayType: Copy,
    T: ArrayContent,
{
    fn as_ref(&self) -> &[u8] {
        self.as_bytes()
    }
}

impl<N, T> Into<String> for ByteArray<N, T>
where
    N: ArrayLength<u8>,
    N::ArrayType: Copy,
    T: ArrayContent,
{
    fn into(self) -> String {
        String::from(self.as_str())
    }
}

impl<N, T> Borrow<[u8]> for ByteArray<N, T>
where
    N: ArrayLength<u8>,
    N::ArrayType: Copy,
    T: ArrayContent,
{
    fn borrow(&self) -> &[u8] {
        self.as_bytes()
    }
}

impl<N, T> Hash for ByteArray<N, T>
where
    N: ArrayLength<u8>,
    N::ArrayType: Copy,
    T: ArrayContent,
{
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.as_bytes().hash(state);
    }
}

impl<N, T> PartialEq for ByteArray<N, T>
where
    N: ArrayLength<u8>,
    N::ArrayType: Copy,
    T: ArrayContent,
{
    fn eq(&self, other: &Self) -> bool {
        self.as_bytes() == other.as_bytes()
    }
}

impl<N, T> Serialize for ByteArray<N, T>
where
    N: ArrayLength<u8>,
    N::ArrayType: Copy,
    T: ArrayContent,
{
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        serializer.serialize_str(self.as_str())
    }
}

impl<'de, N, T> Deserialize<'de> for ByteArray<N, T>
where
    N: ArrayLength<u8>,
    N::ArrayType: Copy,
    T: ArrayContent,
{
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        deserializer.deserialize_str(ByteArrayVisitor {
            phantom_n: PhantomData,
            phantom_t: PhantomData,
        })
    }
}

struct ByteArrayVisitor<N, T> {
    phantom_n: PhantomData<N>,
    phantom_t: PhantomData<T>,
}

impl<'de, N, T> Visitor<'de> for ByteArrayVisitor<N, T>
where
    N: ArrayLength<u8>,
    N::ArrayType: Copy,
    T: ArrayContent,
{
    type Value = ByteArray<N, T>;

    fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        formatter.write_str(T::expected_contents())
    }

    fn visit_str<E>(self, value: &str) -> Result<Self::Value, E>
    where
        E: de::Error,
    {
        Ok(ByteArray::new(value.as_bytes()))
    }
}
