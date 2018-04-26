use fastq::Record;

/// Pointers into a buffer that identify the positions of lines from a FASTQ record
/// header exists at buf[start .. head], seq exists at buf[head .. seq], etc.
#[derive(Deserialize, Serialize, Default, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct ReadOffset {
    exists: bool,
    start: u16,
    head: u16,
    seq: u16,
    qual: u16,
}

/// The possible reads from a Illumina cluster
#[derive(Copy, Clone)]
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

/// Container for all read data from a single Illumina cluster. Faithfully represents
/// the FASTQ data from all available reads, if available.
#[derive(Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Clone)]
pub struct ReadPair {
    offsets: [ReadOffset; 4],

    // Single vector with all the raw FASTQ data
    data: Vec<u8>,
}

impl ReadPair {
    
    pub fn new<R: Record>(rr: [Option<R>; 4]) -> ReadPair {

        let mut offsets = [ReadOffset::default(); 4];
        let mut buf: Vec<u8> = Vec::new();

        for (_rec, which) in rr.iter().zip(WhichRead::read_types().iter()) {
            match _rec {
                &Some(ref rec) => {
                    let start = buf.len() as u16;
                    buf.extend(rec.head());
                    let head = buf.len() as u16;
                    buf.extend(rec.seq());  
                    let seq = buf.len() as u16;
                    buf.extend(rec.head());
                    let qual = buf.len() as u16;
                    let read_offset = ReadOffset { exists: true, start, head, seq, qual };
                    offsets[*which as usize] = read_offset;
                },
                &None => (), // default ReadOffsets is exists = false
            }
        } 

        ReadPair { offsets, data: buf }
    }



    pub fn new_by_part<R: Record>(r1: Option<R>, r2: Option<R>, i1: Option<R>, i2: Option<R>) -> ReadPair {
        Self::new([r1,r2,i1,i2])
    }


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
}

/// Represent a (possibly-correct 10x barcode sequence)
#[derive(Serialize, Deserialize, Clone, Copy, PartialOrd, Ord, PartialEq, Eq, Hash)]
pub enum Barcode {
    Invalid { gem_group: u16, sequence: [u8;16] },
    Valid   { gem_group: u16, sequence: [u8;16] },
}

use self::Barcode::{Valid, Invalid};

impl Barcode {
    pub fn is_valid(&self) -> bool {
        match self {
            &Valid { gem_group, sequence } => true,
            &Invalid { gem_group, sequence } => false,
        }
    }

    pub fn sequence(&self) -> &[u8] {
        match self {
            &Valid { gem_group, ref sequence } => sequence,
            &Invalid { gem_group, ref sequence } => sequence,
        }
    }
}

trait HasBarcode {
    fn barcode(&self) -> Barcode;
}

pub struct Umi {
    length: u8,
    sequence: [u8; 16],
}

trait HasUmi {
    fn umi(&self) -> Umi;
    fn correct_umi(&mut self, corrected_umi: &[u8]);
}


#[derive(Serialize, Deserialize, PartialOrd, Ord, Eq, PartialEq)]
pub struct InlineBcPair {
    data: ReadPair,
    barcode: Barcode,
    read_group_id: u16,
    trim_r1: u8,
    
}

/// A container for the components of a paired-end, barcoded Chromium Genome read.
/// We assume the the 10x barcode is at the beginning of R1
impl InlineBcPair {
    const bc_length: usize = 16;

    pub fn new(data: ReadPair, read_group_id: u16, trim_r1: u8, barcode: Barcode) -> InlineBcPair {
        assert!(data.get(WhichRead::R1, ReadPart::Seq).is_some());
        assert!(data.get(WhichRead::R1, ReadPart::Seq).is_some());

        InlineBcPair {
            data,
            barcode,
            read_group_id,
            trim_r1,
        }
    }

    /// FASTQ read header
    pub fn header(&self) -> &[u8] {
        self.data.get(WhichRead::R1, ReadPart::Header).unwrap()
    }

    /// Full raw R1 sequence
    pub fn r1_seq_raw(&self) -> &[u8] {
        self.data.get(WhichRead::R1, ReadPart::Seq).unwrap()
    }

    /// Full raw R1 QVs
    pub fn r1_qual_raw(&self) -> &[u8] {
        self.data.get(WhichRead::R1, ReadPart::Qual).unwrap()
    }

    /// Full R2 sequence
    pub fn r2_seq(&self) -> &[u8] {
        self.data.get(WhichRead::R2, ReadPart::Seq).unwrap()
    }

    /// Full R2 QVs
    pub fn r2_qual(&self) -> &[u8] {
        self.data.get(WhichRead::R2, ReadPart::Qual).unwrap()
    }

    /// Sample index (I1) sequence
    pub fn si_seq(&self) -> Option<&[u8]> {
        self.data.get(WhichRead::I1, ReadPart::Seq)
    }

    /// Sample index (I1) QVs
    pub fn si_qual(&self) -> Option<&[u8]> {
        self.data.get(WhichRead::I1, ReadPart::Qual)
    }

    /// Raw, uncorrected barcode sequence
    pub fn raw_bc_seq(&self) -> &[u8] {
        &self.r1_seq_raw()[0..Self::bc_length as usize]
    }

    /// Raw barcode QVs
    pub fn raw_bc_qual(&self) -> &[u8] {
        &self.r1_qual_raw()[0..Self::bc_length as usize]
    }

    /// Bases trimmed after the 10x BC, before the start of bases used from R1
    pub fn r1_trim_seq(&self) -> &[u8] {
        &self.r1_seq_raw()[Self::bc_length as usize..(Self::bc_length+self.trim_r1 as usize)]
    }

    /// QVs trimmed after the 10x BC, before the start of bases used from R1
    pub fn r1_trim_qual(&self) -> &[u8] {
        &self.r1_qual_raw()[Self::bc_length as usize..(Self::bc_length+self.trim_r1 as usize)]
    }

    /// Usable R1 bases after removal of BC and trimming
    pub fn r1_seq(&self) -> &[u8] {
        &self.r1_seq_raw()[(Self::bc_length+self.trim_r1 as usize)..]
    }

    /// Usable R1 bases after removal of BC and trimming
    pub fn r1_qual(&self) -> &[u8] {
        &self.r1_qual_raw()[(Self::bc_length+self.trim_r1 as usize)..]
    }
}