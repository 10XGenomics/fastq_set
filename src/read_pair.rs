// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

//! Container for the FASTQ data from a single sequencing 'cluster',
//! including the primary 'R1' and 'R2' and index 'I1' and 'I2' reads.

use bytes::{Bytes, BytesMut};
use failure::Error;
use fastq::{OwnedRecord, Record};
use std::collections::HashMap;
use std::fmt;
use std::io::Write;
use std::ops;
use WhichEnd;

/// Pointers into a buffer that identify the positions of lines from a FASTQ record
/// header exists at buf[start .. head], seq exists at buf[head .. seq], etc.
#[derive(Deserialize, Serialize, Default, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
struct ReadOffset {
    exists: bool,
    start: u16,
    head: u16,
    seq: u16,
    qual: u16,
}

impl ReadOffset {
    fn seq_len(&self) -> Option<usize> {
        if self.exists {
            Some((self.seq - self.head) as usize)
        } else {
            None
        }
    }
}

/// The possible reads from a Illumina cluster. R1 and R2 are the two
/// 'primary' reads, I1 and I2 are the two 'index' samples. I1 is
/// often referred to as the 'sample index read', or I7.  I2 contains
/// the 10x barcode sequence in some 10x assays.
#[derive(Serialize, Deserialize, Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub enum WhichRead {
    R1 = 0,
    R2 = 1,
    I1 = 2,
    I2 = 3,
}

macro_rules! whichread_from {
    ($type:ident) => {
        impl From<$type> for WhichRead {
            fn from(i: $type) -> Self {
                match i {
                    0 => WhichRead::R1,
                    1 => WhichRead::R2,
                    2 => WhichRead::I1,
                    3 => WhichRead::I2,
                    _ => panic!(
                        "Values other than 0,1,2,3 cannot be converted to WhichRead. Got {}",
                        i
                    ),
                }
            }
        }
    };
}
whichread_from!(usize);
whichread_from!(u32);

impl WhichRead {
    pub fn read_types() -> [WhichRead; 4] {
        [WhichRead::R1, WhichRead::R2, WhichRead::I1, WhichRead::I2]
    }
}

/// Components of a FASTQ record.
#[derive(Debug, Copy, Clone)]
pub enum ReadPart {
    Header,
    Seq,
    Qual,
}

/// Compact representation of selected ReadPart and a interval in that read.
/// Supports offsets and lengths up to 32K.
/// Internally it is stored as a `u32` with the following bit layout
///
/// +-----------+--------------+-----------+
/// | WhichRead | Start Offset | Length    |
/// | (2 bits)  | (15 bits)    | (15 bits) |
/// +-----------+--------------+-----------+
/// 
/// Length is optional, with `None` indicating everything until the end of the read.
/// 
/// # Example
/// ```rust
/// extern crate fastq_10x;
/// extern crate fastq;
/// use fastq_10x::read_pair::{RpRange, WhichRead, ReadPart, ReadPair};
/// use fastq_10x::WhichEnd;
/// use fastq::OwnedRecord;
/// let read1 = OwnedRecord {
///     head: b"some_name".to_vec(),
///     seq: b"GTCGCACTGATCTGGGTTAGGCGCGGAGCCGAGGGTTGCACCATTTTTCATTATTGAATGCCAAGATA".to_vec(),
///     qual: b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII".to_vec(),
///     sep: None,
/// };
/// let input = [Some(read1), None, None, None];
/// 
/// // WARNING: DO NOT USE THIS FUNCTION IF YOU ARE STREAMING FASTQ DATA
/// // USE `fastq_10x::read_pair_iter::ReadPairIter` INSTEAD.
/// let read_pair = ReadPair::new(input);
/// // Let's say this read1 is of the form BC(16)-UMI(10)-Insert. This example will
/// // setup different RpRanges to represent these ranges.
/// 
/// let barcode_range = RpRange::new(WhichRead::R1, 0, Some(16));
/// assert_eq!(read_pair.get_range(&barcode_range, ReadPart::Seq).unwrap(), b"GTCGCACTGATCTGGG".to_vec().as_slice());
/// 
/// let umi_range = RpRange::new(WhichRead::R1, 16, Some(10));
/// assert_eq!(read_pair.get_range(&umi_range, ReadPart::Seq).unwrap(), b"TTAGGCGCGG".to_vec().as_slice());
/// 
/// let mut r1_range = RpRange::new(WhichRead::R1, 26, None); // None => everything beyond offset
/// assert_eq!(read_pair.get_range(&r1_range, ReadPart::Seq).unwrap(), b"AGCCGAGGGTTGCACCATTTTTCATTATTGAATGCCAAGATA".to_vec().as_slice());
/// 
/// // Let's say you want to trim first 5 bases in r1
/// r1_range.trim(WhichEnd::FivePrime, 5);
/// assert_eq!(r1_range.offset(), 31);
/// ```
#[derive(Serialize, Deserialize, Copy, Clone, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub struct RpRange {
    val: u32,
}

impl fmt::Debug for RpRange {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "RpRange {{ read: {:?}, offset: {}, len: {:?}, val: {:#b} }}",
            self.read(),
            self.offset(),
            self.len(),
            self.val
        )
    }
}

impl RpRange {
    /// Create a `RpRange` that represent the interval [`offset`, `offset + len`) in
    /// the `read` component of a ReadPair.
    ///
    /// # Args
    /// * `read` - Specify `WhichRead`
    /// * `offset` - Start of the interval. Must be less than 2^15 (=32,768)
    /// * `len` - Optional length that determines the end of the interval. A
    /// value `None` indicates everything from `offset` until the end of the
    /// `read`. Must be less than 2^15 (=32,768)
    ///
    /// # Panics
    /// * If `offset` or `len` is >= `2^15`
    /// 
    /// # Example
    /// ```rust
    /// use fastq_10x::read_pair::RpRange;
    /// use fastq_10x::read_pair::WhichRead;
    /// let range = RpRange::new(WhichRead::R1, 10, Some(50));
    /// assert!(range.read() == WhichRead::R1);
    /// assert!(range.offset() == 10);
    /// assert!(range.len() == Some(50));
    /// ```
    /// 
    /// # Tests
    /// * `test_rprange_invalid_offset()` - Test that this function panics with
    /// an offset that is too large
    /// * `test_rprange_invalid_len()` - Test that this function panics with
    /// a length that is too large
    /// * `prop_test_rprange_representation()` - Test that arbitrary construction of RpRange
    /// stores the values correctly.
    pub fn new(read: WhichRead, offset: usize, len: Option<usize>) -> RpRange {
        assert!(offset < (1 << 15));
        let len_bits = match len {
            Some(v) => {
                assert!(v < (1 << 15));
                v
            }
            None => 0x7FFF,
        };

        let val = (read as u32) << 30 | (offset as u32) << 15 | len_bits as u32;
        RpRange { val }
    }

    #[inline]
    /// Retrieive the read from the internal representation
    pub fn read(&self) -> WhichRead {
        let k = self.val >> 30;
        WhichRead::from(k)
    }

    #[inline]
    /// Retrieive the offset from the internal representation
    pub fn offset(&self) -> usize {
        ((self.val >> 15) & 0x7FFF) as usize
    }

    #[inline]
    /// Retrieive the (optional) length from the internal representation
    pub fn len(&self) -> Option<usize> {
        let len_bits = self.val & 0x7FFF;
        if len_bits == 0x7FFF {
            None
        } else {
            Some(len_bits as usize)
        }
    }

    // Slice the input to the range [offset, offset + len).
    // The input length is not checked in this function.
    fn slice<'a>(&self, input: &'a [u8]) -> &'a [u8] {
        match self.len() {
            Some(l) => &input[self.offset()..self.offset() + l],
            None => &input[self.offset()..],
        }
    }

    // Set the length
    // # Tests
    // * `prop_test_rprange_setter()`
    fn set_len(&mut self, len: usize) {
        assert!(len < (1 << 15));
        self.val = (self.val & (!0x7FFFu32)) | len as u32;
    }

    // Set the offset
    // # Tests
    // * `prop_test_rprange_setter()`
    fn set_offset(&mut self, offset: usize) {
        assert!(offset < (1 << 15));
        self.val = (self.val & (!(0x7FFFu32 << 15))) | (offset as u32) << 15;
    }

    /// Shrink an `RpRange` of known length.
    ///
    /// # Input
    /// * `shrink_range` - `Range<usize>` to shrink the `RpRange` to
    /// 
    /// # Panics
    /// * If the length is set and the `shrink_range` is not within `0..self.len().unwrap()`
    /// * If shrink_range.start > shrink_range.end
    /// 
    /// 
    /// # Example
    /// ```rust
    /// use fastq_10x::read_pair::RpRange;
    /// use fastq_10x::read_pair::WhichRead;
    /// // 100 bases in R1 starting from base 10
    /// let mut rp_range = RpRange::new(WhichRead::R1, 10, Some(100));
    /// let shrink_range = 20..60;
    /// rp_range.shrink(&shrink_range);
    /// assert!(rp_range.read() == WhichRead::R1);
    /// assert!(rp_range.offset() == 30); // 20 + 10
    /// assert!(rp_range.len() == Some(40)); // 60-20
    /// ```
    /// 
    /// # Tests
    /// * `test_shrink_invalid_range_1()`: Test for panic if shrink range start is > length
    /// * `test_shrink_invalid_range_2()`: Test for panic if shrink range end is > length
    /// * `test_shrink_invalid_range_3()`: Test for panic if shrink range start > end
    /// * `test_rprange_trivial_shrink()`: Test shrink to an empty range.
    /// * `prop_test_rprange_shrink()`: Test shrink for arbitrary values of
    /// `RpRange` and valid `shrink_range`
    pub fn shrink(&mut self, shrink_range: &ops::Range<usize>) {

        assert!(
            shrink_range.start <= shrink_range.end,
            "RpRange shrink() expects a valid range with start<=end. Received {:?}", shrink_range
        );

        if let Some(len) = self.len() {
            assert!(
                shrink_range.end <= len,
                "Attempting to shrink more than the current length. shrink_range = {:?}, RpRange = {:?}",
                shrink_range, self,
            );
        }

        let new_offset = self.offset() + shrink_range.start;
        let new_len = shrink_range.end - shrink_range.start;
        self.set_offset(new_offset);
        self.set_len(new_len);
    }

    /// Trim the `RpRange` by specifying the `end` and the `amount`
    /// 
    /// # Inputs
    /// * `end` - whether to trim the `FivePrime` end or the `ThreePrime` end
    /// * `amount` - how many bases to trim.
    /// 
    /// # Panics
    /// * To trim the `ThreePrime` end, the length needs to be known. Panics otherwise
    /// * If length is known, panics if the amount to trim is more than the length
    /// 
    /// # Example
    /// ```rust
    /// use fastq_10x::read_pair::{RpRange, WhichRead};
    /// use fastq_10x::WhichEnd;
    /// let mut rp_range1 = RpRange::new(WhichRead::R1, 40, None);
    /// // Trim 10 bases in the 5' end
    /// rp_range1.trim(WhichEnd::FivePrime, 10);
    /// assert!(rp_range1.read() == WhichRead::R1); // Unchanged
    /// assert!(rp_range1.offset() == 50); // 40 + 10
    /// assert!(rp_range1.len() == None); // Still until the end
    /// 
    /// let mut rp_range2 = RpRange::new(WhichRead::R1, 40, Some(110));
    /// // Trim 20 bases in the 3' end
    /// // Trimming from 3' end needs the length to be set. In this
    /// // case, the range defined is [40, 150)
    /// rp_range2.trim(WhichEnd::ThreePrime, 20);
    /// assert!(rp_range2.read() == WhichRead::R1); // Unchanged
    /// assert!(rp_range2.offset() == 40); // Unchanged
    /// assert!(rp_range2.len() == Some(90)); // 110-20
    /// ```
    /// 
    /// # Tests
    /// * `test_rprange_trim_without_len()` - Make sure 3' trimming without length panics
    /// * `prop_test_rprange_trim()` - Proptest to make sure the trim results are correct
    pub fn trim(&mut self, end: WhichEnd, amount: usize) {

        if amount == 0 { // Trivial case
            return;
        }

        // Panic if we know the length and are asked to trim more than the length
        if let Some(l) = self.len() {
            assert!(
                amount <= l,
                "Attempt to trim more than the length of RpRange"
            );
        }

        match end {
            WhichEnd::ThreePrime => {
                match self.len() {
                    Some(len) => self.set_len(len - amount), // Won't underflow because of the assert above
                    None => panic!("ThreePrime trim is only possible for RpRange with known length!"),
                }
            },
            WhichEnd::FivePrime => {
                let new_offset = self.offset() + amount;
                self.set_offset(new_offset);
                if let Some(l) = self.len() {
                    self.set_len(l - amount); // Won't underflow because of the assert above
                }
            }
        }
    }
}

/// Helper struct used during construction of a ReadPair. The data for the ReadPair is
/// accumulated in the buffer bytes::BytesMut. When all the data has been added, call
/// `freeze()` to convert this into an immutable `ReadPair` object. Multiple `ReadPair` objects
/// will be backed slices into the same buffer, which reduces allocation overhead.
pub(crate) struct MutReadPair<'a> {
    offsets: [ReadOffset; 4],
    data: &'a mut BytesMut,
}

impl<'a> MutReadPair<'a> {
    pub(super) fn empty(buffer: &mut BytesMut) -> MutReadPair {
        let offsets = [ReadOffset::default(); 4];
        MutReadPair {
            offsets,
            data: buffer,
        }
    }

    pub fn new<R: Record>(buffer: &mut BytesMut, rr: [Option<R>; 4]) -> MutReadPair {
        let mut rp = MutReadPair::empty(buffer);

        for (_rec, which) in rr.iter().zip(WhichRead::read_types().iter()) {
            match _rec {
                &Some(ref rec) => rp.push_read(rec, *which),
                &None => (), // default ReadOffsets is exists = false
            }
        }

        rp
    }

    // FIXME: Should we check that the length of seq and qual agree?
    // If we add that check, modify `prop_test_readpair_get()` test
    pub(super) fn push_read<R: Record>(&mut self, rec: &R, which: WhichRead) {
        assert!(!self.offsets[which as usize].exists);
        //let buf = self.data;

        let start = self.data.len() as u16;
        self.data.extend_from_slice(rec.head());
        let head = self.data.len() as u16;
        self.data.extend_from_slice(rec.seq());
        let seq = self.data.len() as u16;
        self.data.extend_from_slice(rec.qual());
        let qual = self.data.len() as u16;
        let read_offset = ReadOffset {
            exists: true,
            start,
            head,
            seq,
            qual,
        };
        self.offsets[which as usize] = read_offset;
    }

    pub fn freeze(self) -> ReadPair {
        ReadPair {
            offsets: self.offsets,
            data: self.data.take().freeze(),
        }
    }
}

/// Container for all read data from a single Illumina cluster. Faithfully represents
/// the FASTQ data from all available reads, if available.
/// Generally should be created by a `ReadPairIter`.
#[derive(Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Clone, Debug)]
pub struct ReadPair {
    offsets: [ReadOffset; 4],

    // Single vector with all the raw FASTQ data.
    // Use with = "serde_bytes" to get much faster perf
    data: Bytes,
}

impl ReadPair {
    #[inline]
    /// Get a ReadPart `part` from a read `which` in this cluster
    pub fn get(&self, which: WhichRead, part: ReadPart) -> Option<&[u8]> {
        if self.offsets[which as usize].exists {
            let w = self.offsets[which as usize];
            match part {
                ReadPart::Header => Some(&self.data[w.start as usize..w.head as usize]),
                ReadPart::Seq => Some(&self.data[w.head as usize..w.seq as usize]),
                ReadPart::Qual => Some(&self.data[w.seq as usize..w.qual as usize]),
            }
        } else {
            None
        }
    }

    #[inline]
    /// Get the range in `RpRange`, return the chosen `part` (sequence or qvs).
    pub fn get_range(&self, rp_range: &RpRange, part: ReadPart) -> Option<&[u8]> {
        let read = self.get(rp_range.read(), part);
        read.map(|r| rp_range.slice(r))
    }

    pub fn to_owned_record(self) -> HashMap<WhichRead, OwnedRecord> {
        let mut result = HashMap::new();
        for &which in WhichRead::read_types().iter() {
            if self.offsets[which as usize].exists {
                let w = self.offsets[which as usize];
                let rec = OwnedRecord {
                    head: self.data[w.start as usize..w.head as usize].to_vec(),
                    seq: self.data[w.head as usize..w.seq as usize].to_vec(),
                    sep: None,
                    qual: self.data[w.seq as usize..w.qual as usize].to_vec(),
                };
                result.insert(which, rec);
            }
        }
        result
    }

    /// Read length of the selected read.
    pub fn len(&self, which: WhichRead) -> Option<usize> {
        self.offsets[which as usize].seq_len()
    }

    /// Write read selected by `which` in FASTQ format to `writer`.
    /// This method will silently do nothing if the selected read doesn't exist.
    pub fn write_fastq<W: Write>(&self, which: WhichRead, writer: &mut W) -> Result<(), Error> {
        if self.offsets[which as usize].exists {
            let head = self.get(which, ReadPart::Header).unwrap();
            writer.write_all(b"@")?;
            writer.write_all(head)?;
            writer.write_all(b"\n")?;

            let seq = self.get(which, ReadPart::Seq).unwrap();
            writer.write_all(seq)?;
            writer.write_all(b"\n+\n")?;

            let qual = self.get(which, ReadPart::Qual).unwrap();
            writer.write_all(qual)?;
            writer.write_all(b"\n")?;
        }

        Ok(())
    }

    /// WARNING: DO NOT USE THIS FUNCTION IF YOU ARE STREAMING FASTQ DATA
    /// This function is intended for testing and illustration purposes
    /// only. Use `ReadPairIter` if you are iterating over a fastq.
    pub fn new<R: Record>(rr: [Option<R>; 4]) -> ReadPair {
        let mut buffer = BytesMut::with_capacity(4096);
        MutReadPair::new(&mut buffer, rr).freeze()
    }
}

#[derive(Debug, Copy, Clone)]
pub struct TrimDef {
    read: WhichRead,
    end: WhichEnd,
    amount: usize,
}

impl TrimDef {
    fn new(read: WhichRead, end: WhichEnd, amount: usize) -> Self {
        TrimDef { read, end, amount }
    }
}

/// A `ReadPair` that supports trimming. All the data corresponding to the
/// original `ReadPair` is retained and the trimming only shifts the
/// offsets for indexing various data parts.
///
/// TODO: Is it better to add trim support to `ReadPair`?
#[derive(Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Clone, Debug)]
pub struct TrimmedReadPair {
    pub trimmed_ranges: [Option<RpRange>; 4],
    readpair: ReadPair,
}

impl From<ReadPair> for TrimmedReadPair {
    fn from(readpair: ReadPair) -> Self {
        let mut trimmed_ranges = [None; 4];
        for &which in WhichRead::read_types().iter() {
            if readpair.offsets[which as usize].exists {
                trimmed_ranges[which as usize] = Some(RpRange::new(
                    which,
                    0usize,
                    readpair.offsets[which as usize].seq_len(),
                ));
            }
        }
        TrimmedReadPair {
            trimmed_ranges,
            readpair,
        }
    }
}

impl TrimmedReadPair {
    /// Trim the reads using the input `TrimDef` which specifies
    /// the read to trim, and how many nucleotides to trim from
    /// which end
    pub fn trim(&mut self, trim_def: TrimDef) -> &mut Self {
        if trim_def.amount == 0 {
            return self;
        }
        assert!(
            self.trimmed_ranges[trim_def.read as usize].is_some(),
            "ERROR: Attempt to trim an empty sequence"
        );
        self.trimmed_ranges[trim_def.read as usize]
            .as_mut()
            .unwrap()
            .trim(trim_def.end, trim_def.amount);
        self
    }

    fn is_trimmed(&self, read: WhichRead) -> bool {
        if let Some(range) = self.trimmed_ranges[read as usize] {
            range.len() == self.readpair.offsets[read as usize].seq_len()
        } else {
            false
        }
    }

    /// Set the length of Illumina read1 to the specified length
    /// If the specified length is greater than the read1 length,
    /// this function does nothing. This function should be called
    /// before applying any other trimming to read1, otherwise, it
    /// will panic.
    pub fn r1_length(&mut self, len: usize) -> &mut Self {
        assert!(
            self.is_trimmed(WhichRead::R1),
            "ERROR: Read1 is already trimmed prior to calling r1_length()"
        );
        assert!(self.readpair.offsets[WhichRead::R1 as usize].exists);
        let ind = WhichRead::R1 as usize;
        let trim_amount = self.readpair.offsets[ind]
            .seq_len()
            .unwrap()
            .saturating_sub(len);
        let trim_def = TrimDef::new(WhichRead::R1, WhichEnd::ThreePrime, trim_amount);
        self.trim(trim_def)
    }

    /// Set the length of Illumina read2 to the specified length
    /// If the specified length is greater than the read2 length,
    /// this function does nothing. This function should be called
    /// before applying any other trimming to read2, otherwise, it
    /// will panic.
    pub fn r2_length(&mut self, len: usize) -> &mut Self {
        assert!(
            self.is_trimmed(WhichRead::R2),
            "ERROR: Read2 is already trimmed prior to calling r2_length()"
        );
        assert!(self.readpair.offsets[WhichRead::R2 as usize].exists);
        let ind = WhichRead::R2 as usize;
        let trim_amount = self.readpair.offsets[ind]
            .seq_len()
            .unwrap()
            .saturating_sub(len);
        let trim_def = TrimDef::new(WhichRead::R2, WhichEnd::ThreePrime, trim_amount);
        self.trim(trim_def)
    }

    #[inline]
    /// Get a ReadPart `part` from a read `which`
    pub fn get(&self, which: WhichRead, part: ReadPart) -> Option<&[u8]> {
        let untrimmed = self.readpair.get(which, part);
        match part {
            ReadPart::Qual | ReadPart::Seq => {
                untrimmed.map(|r| self.trimmed_ranges[which as usize].unwrap().slice(r))
            }
            _ => untrimmed,
        }
    }

    #[inline]
    /// Get the range in `RpRange`, return the chosen `part` (sequence or qvs).
    pub fn get_range(&self, rp_range: &RpRange, part: ReadPart) -> Option<&[u8]> {
        let read = self.get(rp_range.read(), part);
        read.map(|r| rp_range.slice(r))
    }

    pub fn len(&self, which: WhichRead) -> Option<usize> {
        self.trimmed_ranges[which as usize].and_then(|r| r.len())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest;
    use proptest::arbitrary::any;
    use proptest::strategy::Strategy;
    use std::cmp::{max, min};

    const MAX_RPRANGE_ENTRY: usize = 1usize << 15;
    #[test]
    #[should_panic]
    fn test_rprange_invalid_offset() {
        // Offset too large
        let _ = RpRange::new(WhichRead::R1, MAX_RPRANGE_ENTRY, None);
    }

    #[test]
    #[should_panic]
    fn test_rprange_invalid_len() {
        // Length too large
        let _ = RpRange::new(WhichRead::R1, 10, Some(MAX_RPRANGE_ENTRY));
    }

    #[test]
    #[should_panic]
    fn test_rprange_trim_without_len() {
        // 3' Trimming is not allowed if length is unknown
        let mut r = RpRange::new(WhichRead::R1, 10, None);
        r.trim(WhichEnd::ThreePrime, 5);
    }

    #[test]
    fn test_rprange_trivial_shrink() {
        // Trimming is not allowed if length is unknows
        let mut r = RpRange::new(WhichRead::R1, 0, Some(1));
        let shrink_range = 0..0;
        r.shrink(&shrink_range);
        assert_eq!(r.read(), WhichRead::R1);
        assert_eq!(r.offset(), 0);
        assert_eq!(r.len(), Some(0));
    }

    #[test]
    #[should_panic]
    fn test_shrink_invalid_range_1() {
        let mut r = RpRange::new(WhichRead::R1, 10, Some(20));
        r.shrink(&(40..50));
    }

    #[test]
    #[should_panic]
    fn test_shrink_invalid_range_2() {
        let mut r = RpRange::new(WhichRead::R1, 10, Some(20));
        r.shrink(&(10..50));
    }

    #[test]
    #[should_panic]
    fn test_shrink_invalid_range_3() {
        let mut r = RpRange::new(WhichRead::R1, 10, Some(20));
        r.shrink(&(10..5));
    }

    proptest! {
        #[test]
        fn prop_test_rprange_representation(
            read in (0..4usize).prop_map(|r| WhichRead::from(r)),
            offset in 0..MAX_RPRANGE_ENTRY,
            len in any::<usize>()
        ) {
            // Make sure out internal compact representation is valid for random inputs
            let len = if len <= MAX_RPRANGE_ENTRY { Some(len) } else { None };
            let rprange = RpRange::new(read, offset, len);
            assert!(rprange.read() == read);
            assert!(rprange.offset() == offset);
            assert!(rprange.len() == len);
        }
    }

    proptest! {
        #[test]
        fn prop_test_rprange_setter(
            read in (0..4usize).prop_map(|r| WhichRead::from(r)),
            offset in 0..MAX_RPRANGE_ENTRY,
            len in any::<usize>(),
            new_offset in 0..MAX_RPRANGE_ENTRY,
            new_len in 0..MAX_RPRANGE_ENTRY
        ) {
            // Make sure we set length and offset correctly using the setter functions
            // for random inputs
            let len = if len <= MAX_RPRANGE_ENTRY { Some(len) } else { None };
            let mut rprange = RpRange::new(read, offset, len);
            rprange.set_offset(new_offset);
            rprange.set_len(new_len);
            assert!(rprange.read() == read);
            assert!(rprange.offset() == new_offset);
            assert!(rprange.len() == Some(new_len));
        }
    }

    proptest! {
        #[test]
        fn prop_test_rprange_trim(
            val in any::<u32>(),
            end in any::<bool>(),
            amount in any::<usize>()
        ) {
            // Make sure we trim the range correctly
            let mut rprange = RpRange { val };
            let old_offset = rprange.offset();
            let old_read = rprange.read();
            if let Some(l) = rprange.len() {
                let amount = amount % max(l, 1);
                if amount + old_offset < MAX_RPRANGE_ENTRY {
                    if end {
                        rprange.trim(WhichEnd::ThreePrime, amount);
                        assert_eq!(rprange.read(), old_read);
                        assert_eq!(rprange.offset(), old_offset);
                        assert_eq!(rprange.len(), Some(l-amount));
                    } else {
                        rprange.trim(WhichEnd::FivePrime, amount);
                        assert_eq!(rprange.read(), old_read);
                        assert_eq!(rprange.offset(), old_offset + amount);
                        assert_eq!(rprange.len(), Some(l-amount));
                    }
                }
            } else {
                rprange.trim(WhichEnd::FivePrime, amount);
                assert_eq!(rprange.read(), old_read);
                assert_eq!(rprange.offset(), old_offset + amount);
                assert_eq!(rprange.len(), None);
            }

        }
    }

    proptest! {
        #[test]
        fn prop_test_rprange_shrink(
            val in any::<u32>(),
            x in any::<usize>(),
            y in any::<usize>()
        ) {
            // Make sure we shrink the range correctly
            let mut rprange = RpRange { val };

            let max_val = rprange.len().unwrap_or(MAX_RPRANGE_ENTRY);
            let x = min(x, max_val);
            let y = min(y, max_val);
            let shrink_range = min(x, y)..max(x, y);

            let old_offset = rprange.offset();
            let old_read = rprange.read();

            let expected_offset = old_offset + shrink_range.start;
            let expected_len = Some(shrink_range.end - shrink_range.start);
            if expected_offset < MAX_RPRANGE_ENTRY {
                rprange.shrink(&shrink_range);
                assert_eq!(rprange.read(), old_read);
                assert_eq!(rprange.offset(), expected_offset);
                assert_eq!(rprange.len(), expected_len);
            }

        }
    }

    proptest! {
        #[test]
        fn prop_test_readpair_get(
            ref head in proptest::collection::vec(any::<u8>(), 0usize..500usize),
            ref seq in proptest::collection::vec(any::<u8>(), 0usize..1000usize),
            ref qual in proptest::collection::vec(any::<u8>(), 0usize..1000usize),
            pos in 0..4usize
        ) {
            let mut buffer = BytesMut::with_capacity(4096);
            let owned = OwnedRecord {
                head: head.clone(),
                seq: seq.clone(),
                qual: qual.clone(),
                sep: None,
            };
            let mut input = [None, None, None, None];
            input[pos] = Some(owned);
            let read_pair = MutReadPair::new(&mut buffer, input).freeze();
            let read = WhichRead::from(pos);
            assert_eq!(read_pair.get(read, ReadPart::Header), Some(head.as_slice()));
            assert_eq!(read_pair.get(read, ReadPart::Qual), Some(qual.as_slice()));
            assert_eq!(read_pair.get(read, ReadPart::Seq), Some(seq.as_slice()));
        }
    }
}
