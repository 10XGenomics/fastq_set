// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

use fastq::Record;

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

pub enum ReadPart {
    Header,
    Seq,
    Qual,
}

/// A compact representation of a slice into a ReadPair.
/// Use this when you need to store a 
#[derive(Serialize, Deserialize, Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub struct RpRange {
    val: u32
}

impl RpRange {
    pub fn new(read: WhichRead, offset: usize, len: Option<usize>) -> RpRange {
        assert!(offset < (1<<15));
        let len_bits = match len {
            Some(v) => {
                assert!(v < (1<<15));
                v
            },
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
            _ => unreachable!("bad read id")
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

    pub fn slice<'a>(&self, input: &'a[u8]) -> &'a[u8] {
        match self.len() {
            Some(l) => &input[self.offset() .. self.offset() + l],
            None => &input[self.offset() .. ],
        }
    }
}


/// Container for all read data from a single Illumina cluster. Faithfully represents
/// the FASTQ data from all available reads, if available.
#[derive(Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Clone, Debug)]
pub struct ReadPair {
    offsets: [ReadOffset; 4],

    // Single vector with all the raw FASTQ data
    data: Vec<u8>,
}

impl ReadPair {

    pub(super) fn empty() -> ReadPair {
        let offsets = [ReadOffset::default(); 4];
        let data = Vec::new();
        ReadPair { offsets, data }
    }

    fn new<R: Record>(rr: [Option<R>; 4]) -> ReadPair {
        let offsets = [ReadOffset::default(); 4];
        let data = Vec::new();
        let mut rp = ReadPair { offsets, data };

        for (_rec, which) in rr.iter().zip(WhichRead::read_types().iter()) {
            match _rec {
                &Some(ref rec) => {
                    rp.push_read(rec, *which)
                },
                &None => (), // default ReadOffsets is exists = false
            }
        } 

        rp
    }

    pub fn new_by_part<R: Record>(r1: Option<R>, r2: Option<R>, i1: Option<R>, i2: Option<R>) -> ReadPair {
        Self::new([r1,r2,i1,i2])
    }

    pub(super) fn push_read<R:Record>(&mut self, rec: &R, which: WhichRead) {
        let buf = &mut self.data;

        let start = buf.len() as u16;
        buf.extend(rec.head());
        let head = buf.len() as u16;
        buf.extend(rec.seq());  
        let seq = buf.len() as u16;
        buf.extend(rec.head());
        let qual = buf.len() as u16;
        let read_offset = ReadOffset { exists: true, start, head, seq, qual };
        self.offsets[which as usize] = read_offset;
    }


    #[inline]
    /// Get a ReadPart `part` from a read `which` in this cluster
    pub fn get(&self, which: WhichRead, part: ReadPart) -> Option<&[u8]> {
        if self.offsets[which as usize].exists {
            let w = self.offsets[which as usize];
            match part {
                ReadPart::Header => Some(&self.data[w.start as usize .. w.head as usize]),
                ReadPart::Seq => Some(&self.data[w.head as usize .. w.seq as usize]),
                ReadPart::Qual => Some(&self.data[w.seq as usize .. w.qual as usize]),
            }
        } else {
            None
        }
    }

    #[inline]
    pub fn get_range(&self, rp_range: &RpRange, part: ReadPart) -> Option<&[u8]> {
        let read = self.get(rp_range.read(), part);
        read.map(|r| rp_range.slice(r))
    }
}