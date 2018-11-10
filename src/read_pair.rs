// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

//! Container for the FASTQ data from a single sequencing 'cluster',
//! including the primary 'R1' and 'R2' and index 'I1' and 'I2' reads.

use fastq::{Record, OwnedRecord};
use std::collections::HashMap;
use std::fmt;
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

#[derive(Serialize, Deserialize, Copy, Clone, Eq, PartialEq, PartialOrd, Ord, Hash)]
pub struct RpRange {
    val: u32,
}

impl fmt::Debug for RpRange {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "RpRange {{ read: {:?}, offset: {}, len: {:?}, val: {:#b} }}", self.read(), self.offset(), self.len(), self.val)
    }
}

impl RpRange {
    /// Create a `RpRange` that represent the interval [offset, offset_length) in
    /// the `read` component of a ReadPair.
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
    pub fn read(&self) -> WhichRead {
        let k = self.val >> 30;
        match k {
            0 => WhichRead::R1,
            1 => WhichRead::R2,
            2 => WhichRead::I1,
            3 => WhichRead::I2,
            _ => unreachable!("bad read id"),
        }
    }

    #[inline]
    pub fn offset(&self) -> usize {
        ((self.val >> 15) & 0x7FFF) as usize
    }

    #[inline]
    pub fn len(&self) -> Option<usize> {
        let len_bits = self.val & 0x7FFF;
        if len_bits == 0x7FFF {
            None
        } else {
            Some(len_bits as usize)
        }
    }

    fn slice<'a>(&self, input: &'a [u8]) -> &'a [u8] {
        match self.len() {
            Some(l) => &input[self.offset()..self.offset() + l],
            None => &input[self.offset()..],
        }
    }

    pub fn set_len(&mut self, len: usize) {
        assert!(len < (1 << 15));
        self.val = (self.val & (!0x7FFFu32)) | len as u32;
    }

    fn set_offset(&mut self, offset: usize) {
        assert!(offset < (1 << 15));
        self.val = (self.val & (!(0x7FFFu32 << 15))) | (offset as u32) << 15;
    }

    pub fn trim(&mut self, trim_def: TrimDef) {
        assert!(trim_def.read==self.read());

        if trim_def.amount > 0 {
            match trim_def.end {
                WhichEnd::ThreePrime => {
                    if let Some(l) = self.len() {
                        self.set_len(l.saturating_sub(trim_def.amount));
                    } else {
                        panic!("Cannot apply trim to the 3' end without knowing the length");
                    }
                },
                WhichEnd::FivePrime => {
                    let offset = self.offset();
                    self.set_offset(offset + trim_def.amount);
                    if let Some(l) = self.len() {
                        self.set_len(l.saturating_sub(trim_def.amount));
                    }
                }
            }
        }
    }
}

use serde_bytes::ByteBuf;

/// Container for all read data from a single Illumina cluster. Faithfully represents
/// the FASTQ data from all available reads, if available.
/// Generally should be created by a `ReadPairIter`.
#[derive(Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Clone, Debug)]
pub struct ReadPair {
    offsets: [ReadOffset; 4],

    // Single vector with all the raw FASTQ data
    // Use with = "serde_bytes" to get much faster perf
    data: ByteBuf,
}

impl ReadPair {
    // Make space for the full read pair in one allocation
    const RP_CAPACITY: usize = 1024;

    pub(super) fn empty() -> ReadPair {
        let offsets = [ReadOffset::default(); 4];
        let data = Vec::with_capacity(Self::RP_CAPACITY);
        ReadPair {
            offsets,
            data: ByteBuf::from(data),
        }
    }

    pub fn new<R: Record>(rr: [Option<R>; 4]) -> ReadPair {
        let offsets = [ReadOffset::default(); 4];
        let data = Vec::with_capacity(Self::RP_CAPACITY);
        let mut rp = ReadPair {
            offsets,
            data: ByteBuf::from(data),
        };

        for (_rec, which) in rr.iter().zip(WhichRead::read_types().iter()) {
            match _rec {
                &Some(ref rec) => rp.push_read(rec, *which),
                &None => (), // default ReadOffsets is exists = false
            }
        }

        rp
    }

    pub(super) fn push_read<R: Record>(&mut self, rec: &R, which: WhichRead) {
        let buf: &mut Vec<u8> = self.data.as_mut();

        let start = buf.len() as u16;
        buf.extend(rec.head());
        let head = buf.len() as u16;
        buf.extend(rec.seq());
        let seq = buf.len() as u16;
        buf.extend(rec.qual());
        let qual = buf.len() as u16;
        let read_offset = ReadOffset {
            exists: true,
            start,
            head,
            seq,
            qual,
        };
        self.offsets[which as usize] = read_offset;
    }

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

    pub fn len(&self, which: WhichRead) -> Option<usize> {
        self.offsets[which as usize].seq_len()
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
        TrimDef {
            read,
            end,
            amount,
        }
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
                trimmed_ranges[which as usize] = Some(RpRange::new(which, 0usize, readpair.offsets[which as usize].seq_len()));
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
        assert!(self.trimmed_ranges[trim_def.read as usize].is_some(), "ERROR: Attempt to trim an empty sequence");
        self.trimmed_ranges[trim_def.read as usize].as_mut().unwrap().trim(trim_def);
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
        assert!(self.is_trimmed(WhichRead::R1), "ERROR: Read1 is already trimmed prior to calling r1_length()");
        assert!(self.readpair.offsets[WhichRead::R1 as usize].exists);
        let ind = WhichRead::R1 as usize;
        let trim_amount = self.readpair.offsets[ind].seq_len().unwrap().saturating_sub(len);
        let trim_def = TrimDef::new(WhichRead::R1, WhichEnd::ThreePrime, trim_amount);
        self.trim(trim_def)
    }

    /// Set the length of Illumina read2 to the specified length
    /// If the specified length is greater than the read2 length,
    /// this function does nothing. This function should be called
    /// before applying any other trimming to read2, otherwise, it
    /// will panic.
    pub fn r2_length(&mut self, len: usize) -> &mut Self {
        assert!(self.is_trimmed(WhichRead::R2), "ERROR: Read2 is already trimmed prior to calling r2_length()");
        assert!(self.readpair.offsets[WhichRead::R2 as usize].exists);
        let ind = WhichRead::R2 as usize;
        let trim_amount = self.readpair.offsets[ind].seq_len().unwrap().saturating_sub(len);
        let trim_def = TrimDef::new(WhichRead::R2, WhichEnd::ThreePrime, trim_amount);
        self.trim(trim_def)
    }

    #[inline]
    /// Get a ReadPart `part` from a read `which`
    pub fn get(&self, which: WhichRead, part: ReadPart) -> Option<&[u8]> {
        let untrimmed = self.readpair.get(which, part);
        match part {
            ReadPart::Qual | ReadPart::Seq => untrimmed.map(|r| self.trimmed_ranges[which as usize].unwrap().slice(r)),
            _ => untrimmed
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