
use {HasBarcode, Barcode, FastqProcessor, AlignableRead, SampleDef, FastqFiles, SSeq};
use raw::{ReadPair, WhichRead, ReadPart, RpRange};

#[derive(Serialize, Deserialize, PartialOrd, PartialEq, Debug, Clone)]
pub struct DnaChunk {
    barcode: Option<String>,
    barcode_reverse_complement: bool,
    bc_in_read: Option<u8>,
    bc_length: Option<usize>,
    gem_group: u16,
    read1: String,
    read2: Option<String>,
    read_group: String,
    reads_interleaved: bool,
    sample_index: Option<String>,
    subsample_rate: f64,
    read_group_id: Option<u16>
}

#[derive(Serialize, Deserialize, PartialOrd, Ord, PartialEq, Eq)]
pub struct DnaFastqChunk {
    bc_in_read: Option<usize>,
    gem_group: usize,
    read_group: String,
    read_group_id: u16,
    trim_r1: u8,
    sample_def: SampleDef,
}

impl FastqProcessor<DnaRead> for DnaChunk {
    fn process_read(&self, read: ReadPair) -> Option<DnaRead> {
        assert!(read.get(WhichRead::R1, ReadPart::Seq).is_some());
        assert!(read.get(WhichRead::R2, ReadPart::Seq).is_some());

        // Setup initial (uncorrected) bacode
        let bc_range = match self.bc_in_read {
            Some(r) => RpRange::new(WhichRead::R1, 0, Some(r as usize)),
            None => RpRange::new(WhichRead::I2, 0, None),
        };
        let barcode = Barcode::new(self.gem_group, read.get_range(&bc_range, ReadPart::Seq).unwrap(), true);

        Some(DnaRead {
            data: read,
            barcode,
            read_group_id: self.read_group_id.expect("need to set read group id"),
            trim_r1: 8
        })
    }

    fn description(&self) -> String {
        self.read1.clone()
    }

    fn fastq_files(&self) -> FastqFiles {
        FastqFiles {
            r1: self.read1.clone(),
            r2: self.read2.clone(),
            i1: self.sample_index.clone(),
            i2: self.barcode.clone(),
            r1_interleaved: self.reads_interleaved,
        }
    }

    fn bc_subsample_rate(&self) -> f64 {
        1.0
    }

    fn read_subsample_rate(&self) -> f64 {
        1.0
    }
}

#[derive(Serialize, Deserialize, PartialOrd, Ord, Eq, PartialEq)]
pub struct DnaRead {
    data: ReadPair,
    barcode: Barcode,
    trim_r1: u8,
    read_group_id: u16,
}

impl HasBarcode for DnaRead {
    fn barcode(&self) -> Barcode {
        self.barcode
    }
}

/// A container for the components of a paired-end, barcoded Chromium Genome read.
/// We assume the the 10x barcode is at the beginning of R1
impl DnaRead {
    const bc_length: usize = 16;


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


impl AlignableRead for DnaRead {
    fn alignable_reads(&self) -> (Option<&[u8]>, Option<&[u8]>) {
        (Some(self.r1_seq()), Some(self.r2_seq()))
    }
}


#[cfg(test)]
mod test_dna_cfg {
    use super::*;
    use serde_json;

    fn load_dna_chunk_def(chunk_json: &str) -> Vec<DnaChunk> {
        let mut c: Vec<DnaChunk> = serde_json::from_str(chunk_json).unwrap();

        for i in 0 .. c.len() {
            c[i].read_group_id = Some(i as u16);
        }

        c
    }

    #[test]
    fn test_crg_cfg() {
        let c = load_dna_chunk_def(CRG_CFG);
        println!("{:?}", c);
    }

    #[test]
    fn test_atac_cfg() {
        let c = load_dna_chunk_def(ATAC_CFG);
        println!("{:?}", c);
    }

    #[test]
    fn test_scdna_cfg() {
        let c = load_dna_chunk_def(SCDNA_CFG);
        println!("{:?}", c);
    }

    #[test]
    fn test_load_atac() {
        let c = load_dna_chunk_def(ATAC_CFG_TEST);
        println!("{:?}", c);
        let sorter = ::bc_sort::SortAndCorrect::<_,DnaRead>::new(c);
        sorter.sort_bcs("test/rp_atac/reads", "test/rp_atac/bc_counts", 2).unwrap();
    }

    const CRG_CFG: &str = r#"
    [
        {
            "barcode": null,
            "barcode_reverse_complement": false,
            "bc_in_read": 1,
            "bc_length": 16,
            "gem_group": 1,
            "read1": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-RA_si-AAGGGTGA_lane-001-chunk-002.fastq.gz",
            "read2": null,
            "read_group": "20486:MissingLibrary:1:unknown_fc:0",
            "reads_interleaved": true,
            "sample_index": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-I1_si-AAGGGTGA_lane-001-chunk-002.fastq.gz",
            "subsample_rate": 0.3224737982551789
        },
        {
            "barcode": null,
            "barcode_reverse_complement": false,
            "bc_in_read": 1,
            "bc_length": 16,
            "gem_group": 1,
            "read1": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-RA_si-AAGGGTGA_lane-002-chunk-004.fastq.gz",
            "read2": null,
            "read_group": "20486:MissingLibrary:1:unknown_fc:0",
            "reads_interleaved": true,
            "sample_index": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-I1_si-AAGGGTGA_lane-002-chunk-004.fastq.gz",
            "subsample_rate": 0.3224737982551789
        },
        {
            "barcode": null,
            "barcode_reverse_complement": false,
            "bc_in_read": 1,
            "bc_length": 16,
            "gem_group": 1,
            "read1": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-RA_si-AAGGGTGA_lane-003-chunk-007.fastq.gz",
            "read2": null,
            "read_group": "20486:MissingLibrary:1:unknown_fc:0",
            "reads_interleaved": true,
            "sample_index": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-I1_si-AAGGGTGA_lane-003-chunk-007.fastq.gz",
            "subsample_rate": 0.3224737982551789
        },
        {
            "barcode": null,
            "barcode_reverse_complement": false,
            "bc_in_read": 1,
            "bc_length": 16,
            "gem_group": 1,
            "read1": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-RA_si-AAGGGTGA_lane-004-chunk-006.fastq.gz",
            "read2": null,
            "read_group": "20486:MissingLibrary:1:unknown_fc:0",
            "reads_interleaved": true,
            "sample_index": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-I1_si-AAGGGTGA_lane-004-chunk-006.fastq.gz",
            "subsample_rate": 0.3224737982551789
        },
        {
            "barcode": null,
            "barcode_reverse_complement": false,
            "bc_in_read": 1,
            "bc_length": 16,
            "gem_group": 1,
            "read1": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-RA_si-AAGGGTGA_lane-005-chunk-000.fastq.gz",
            "read2": null,
            "read_group": "20486:MissingLibrary:1:unknown_fc:0",
            "reads_interleaved": true,
            "sample_index": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-I1_si-AAGGGTGA_lane-005-chunk-000.fastq.gz",
            "subsample_rate": 0.3224737982551789
        },
        {
            "barcode": null,
            "barcode_reverse_complement": false,
            "bc_in_read": 1,
            "bc_length": 16,
            "gem_group": 1,
            "read1": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-RA_si-AAGGGTGA_lane-006-chunk-003.fastq.gz",
            "read2": null,
            "read_group": "20486:MissingLibrary:1:unknown_fc:0",
            "reads_interleaved": true,
            "sample_index": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-I1_si-AAGGGTGA_lane-006-chunk-003.fastq.gz",
            "subsample_rate": 0.3224737982551789
        },
        {
            "barcode": null,
            "barcode_reverse_complement": false,
            "bc_in_read": 1,
            "bc_length": 16,
            "gem_group": 1,
            "read1": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-RA_si-AAGGGTGA_lane-007-chunk-005.fastq.gz",
            "read2": null,
            "read_group": "20486:MissingLibrary:1:unknown_fc:0",
            "reads_interleaved": true,
            "sample_index": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-I1_si-AAGGGTGA_lane-007-chunk-005.fastq.gz",
            "subsample_rate": 0.3224737982551789
        },
        {
            "barcode": null,
            "barcode_reverse_complement": false,
            "bc_in_read": 1,
            "bc_length": 16,
            "gem_group": 1,
            "read1": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-RA_si-AAGGGTGA_lane-008-chunk-001.fastq.gz",
            "read2": null,
            "read_group": "20486:MissingLibrary:1:unknown_fc:0",
            "reads_interleaved": true,
            "sample_index": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-I1_si-AAGGGTGA_lane-008-chunk-001.fastq.gz",
            "subsample_rate": 0.3224737982551789
        },
        {
            "barcode": null,
            "barcode_reverse_complement": false,
            "bc_in_read": 1,
            "bc_length": 16,
            "gem_group": 1,
            "read1": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-RA_si-NAGGGTGA_lane-001-chunk-002.fastq.gz",
            "read2": null,
            "read_group": "20486:MissingLibrary:1:unknown_fc:0",
            "reads_interleaved": true,
            "sample_index": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-I1_si-NAGGGTGA_lane-001-chunk-002.fastq.gz",
            "subsample_rate": 0.3224737982551789
        },
        {
            "barcode": null,
            "barcode_reverse_complement": false,
            "bc_in_read": 1,
            "bc_length": 16,
            "gem_group": 1,
            "read1": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-RA_si-NAGGGTGA_lane-002-chunk-004.fastq.gz",
            "read2": null,
            "read_group": "20486:MissingLibrary:1:unknown_fc:0",
            "reads_interleaved": true,
            "sample_index": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-I1_si-NAGGGTGA_lane-002-chunk-004.fastq.gz",
            "subsample_rate": 0.3224737982551789
        },
        {
            "barcode": null,
            "barcode_reverse_complement": false,
            "bc_in_read": 1,
            "bc_length": 16,
            "gem_group": 1,
            "read1": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-RA_si-NAGGGTGA_lane-003-chunk-007.fastq.gz",
            "read2": null,
            "read_group": "20486:MissingLibrary:1:unknown_fc:0",
            "reads_interleaved": true,
            "sample_index": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-I1_si-NAGGGTGA_lane-003-chunk-007.fastq.gz",
            "subsample_rate": 0.3224737982551789
        },
        {
            "barcode": null,
            "barcode_reverse_complement": false,
            "bc_in_read": 1,
            "bc_length": 16,
            "gem_group": 1,
            "read1": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-RA_si-NAGGGTGA_lane-004-chunk-006.fastq.gz",
            "read2": null,
            "read_group": "20486:MissingLibrary:1:unknown_fc:0",
            "reads_interleaved": true,
            "sample_index": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-I1_si-NAGGGTGA_lane-004-chunk-006.fastq.gz",
            "subsample_rate": 0.3224737982551789
        },
        {
            "barcode": null,
            "barcode_reverse_complement": false,
            "bc_in_read": 1,
            "bc_length": 16,
            "gem_group": 1,
            "read1": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-RA_si-NAGGGTGA_lane-005-chunk-000.fastq.gz",
            "read2": null,
            "read_group": "20486:MissingLibrary:1:unknown_fc:0",
            "reads_interleaved": true,
            "sample_index": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-I1_si-NAGGGTGA_lane-005-chunk-000.fastq.gz",
            "subsample_rate": 0.3224737982551789
        },
        {
            "barcode": null,
            "barcode_reverse_complement": false,
            "bc_in_read": 1,
            "bc_length": 16,
            "gem_group": 1,
            "read1": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-RA_si-NAGGGTGA_lane-006-chunk-003.fastq.gz",
            "read2": null,
            "read_group": "20486:MissingLibrary:1:unknown_fc:0",
            "reads_interleaved": true,
            "sample_index": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-I1_si-NAGGGTGA_lane-006-chunk-003.fastq.gz",
            "subsample_rate": 0.3224737982551789
        },
        {
            "barcode": null,
            "barcode_reverse_complement": false,
            "bc_in_read": 1,
            "bc_length": 16,
            "gem_group": 1,
            "read1": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-RA_si-NAGGGTGA_lane-007-chunk-005.fastq.gz",
            "read2": null,
            "read_group": "20486:MissingLibrary:1:unknown_fc:0",
            "reads_interleaved": true,
            "sample_index": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-I1_si-NAGGGTGA_lane-007-chunk-005.fastq.gz",
            "subsample_rate": 0.3224737982551789
        },
        {
            "barcode": null,
            "barcode_reverse_complement": false,
            "bc_in_read": 1,
            "bc_length": 16,
            "gem_group": 1,
            "read1": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-RA_si-NAGGGTGA_lane-008-chunk-001.fastq.gz",
            "read2": null,
            "read_group": "20486:MissingLibrary:1:unknown_fc:0",
            "reads_interleaved": true,
            "sample_index": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-I1_si-NAGGGTGA_lane-008-chunk-001.fastq.gz",
            "subsample_rate": 0.3224737982551789
        },
        {
            "barcode": null,
            "barcode_reverse_complement": false,
            "bc_in_read": 1,
            "bc_length": 16,
            "gem_group": 1,
            "read1": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-RA_si-NTACTCAG_lane-001-chunk-002.fastq.gz",
            "read2": null,
            "read_group": "20486:MissingLibrary:1:unknown_fc:0",
            "reads_interleaved": true,
            "sample_index": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-I1_si-NTACTCAG_lane-001-chunk-002.fastq.gz",
            "subsample_rate": 0.3224737982551789
        },
        {
            "barcode": null,
            "barcode_reverse_complement": false,
            "bc_in_read": 1,
            "bc_length": 16,
            "gem_group": 1,
            "read1": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-RA_si-NTACTCAG_lane-002-chunk-004.fastq.gz",
            "read2": null,
            "read_group": "20486:MissingLibrary:1:unknown_fc:0",
            "reads_interleaved": true,
            "sample_index": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-I1_si-NTACTCAG_lane-002-chunk-004.fastq.gz",
            "subsample_rate": 0.3224737982551789
        },
        {
            "barcode": null,
            "barcode_reverse_complement": false,
            "bc_in_read": 1,
            "bc_length": 16,
            "gem_group": 1,
            "read1": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-RA_si-NTACTCAG_lane-003-chunk-007.fastq.gz",
            "read2": null,
            "read_group": "20486:MissingLibrary:1:unknown_fc:0",
            "reads_interleaved": true,
            "sample_index": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-I1_si-NTACTCAG_lane-003-chunk-007.fastq.gz",
            "subsample_rate": 0.3224737982551789
        },
        {
            "barcode": null,
            "barcode_reverse_complement": false,
            "bc_in_read": 1,
            "bc_length": 16,
            "gem_group": 1,
            "read1": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-RA_si-NTACTCAG_lane-004-chunk-006.fastq.gz",
            "read2": null,
            "read_group": "20486:MissingLibrary:1:unknown_fc:0",
            "reads_interleaved": true,
            "sample_index": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-I1_si-NTACTCAG_lane-004-chunk-006.fastq.gz",
            "subsample_rate": 0.3224737982551789
        },
        {
            "barcode": null,
            "barcode_reverse_complement": false,
            "bc_in_read": 1,
            "bc_length": 16,
            "gem_group": 1,
            "read1": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-RA_si-NTACTCAG_lane-005-chunk-000.fastq.gz",
            "read2": null,
            "read_group": "20486:MissingLibrary:1:unknown_fc:0",
            "reads_interleaved": true,
            "sample_index": "/mnt/analysis/marsoc/pipestances/HMVT3CCXX/BCL_PROCESSOR_PD/HMVT3CCXX/1015.9.4-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/DEMULTIPLEX/fork0/files/demultiplexed_fastq_path/read-I1_si-NTACTCAG_lane-005-chunk-000.fastq.gz",
            "subsample_rate": 0.3224737982551789
        }
    ]"#;

    const ATAC_CFG_TEST: &str = r#"
    [
        {
            "barcode": "test/bcl_processor/atac/read-I2_si-CGGAGCAC_lane-001-chunk-001.fastq.gz",
            "barcode_reverse_complement": false,
            "gem_group": 1,
            "read1": "test/bcl_processor/atac/read-RA_si-CGGAGCAC_lane-001-chunk-001.fastq.gz",
            "read2": null,
            "read_group": "66333:66333:1:HC7WVDMXX:1",
            "reads_interleaved": true,
            "sample_index": "test/bcl_processor/atac/read-I1_si-CGGAGCAC_lane-001-chunk-001.fastq.gz",
            "subsample_rate": 0.37474273634185645
        },
        {
            "barcode": "test/bcl_processor/atac/read-I2_si-GACCTATT_lane-001-chunk-001.fastq.gz",
            "barcode_reverse_complement": false,
            "gem_group": 1,
            "read1": "test/bcl_processor/atac/read-RA_si-GACCTATT_lane-001-chunk-001.fastq.gz",
            "read2": null,
            "read_group": "66333:66333:1:HC7WVDMXX:1",
            "reads_interleaved": true,
            "sample_index": "test/bcl_processor/atac/read-I1_si-GACCTATT_lane-001-chunk-001.fastq.gz",
            "subsample_rate": 0.37474273634185645
        }
    ]
    "#;


    const ATAC_CFG: &str = r#"
    [
        {
            "barcode": "/mnt/analysis/marsoc/pipestances/HC7WVDMXX/BCL_PROCESSOR_PD/HC7WVDMXX/1015.12.13-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/MERGE_FASTQS_FROM_TILES/fork0/join-u37c6ba9993/files/demultiplexed_fastq_path/read-I2_si-CGGAGCAC_lane-001-chunk-001.fastq.gz",
            "barcode_reverse_complement": false,
            "gem_group": 1,
            "read1": "/mnt/analysis/marsoc/pipestances/HC7WVDMXX/BCL_PROCESSOR_PD/HC7WVDMXX/1015.12.13-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/MERGE_FASTQS_FROM_TILES/fork0/join-u37c6ba9993/files/demultiplexed_fastq_path/read-RA_si-CGGAGCAC_lane-001-chunk-001.fastq.gz",
            "read2": null,
            "read_group": "66333:66333:1:HC7WVDMXX:1",
            "reads_interleaved": true,
            "sample_index": "/mnt/analysis/marsoc/pipestances/HC7WVDMXX/BCL_PROCESSOR_PD/HC7WVDMXX/1015.12.13-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/MERGE_FASTQS_FROM_TILES/fork0/join-u37c6ba9993/files/demultiplexed_fastq_path/read-I1_si-CGGAGCAC_lane-001-chunk-001.fastq.gz",
            "subsample_rate": 0.37474273634185645
        },
        {
            "barcode": "/mnt/analysis/marsoc/pipestances/HC7WVDMXX/BCL_PROCESSOR_PD/HC7WVDMXX/1015.12.13-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/MERGE_FASTQS_FROM_TILES/fork0/join-u37c6ba9993/files/demultiplexed_fastq_path/read-I2_si-CGGAGCAC_lane-002-chunk-001.fastq.gz",
            "barcode_reverse_complement": false,
            "gem_group": 1,
            "read1": "/mnt/analysis/marsoc/pipestances/HC7WVDMXX/BCL_PROCESSOR_PD/HC7WVDMXX/1015.12.13-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/MERGE_FASTQS_FROM_TILES/fork0/join-u37c6ba9993/files/demultiplexed_fastq_path/read-RA_si-CGGAGCAC_lane-002-chunk-001.fastq.gz",
            "read2": null,
            "read_group": "66333:66333:1:HC7WVDMXX:2",
            "reads_interleaved": true,
            "sample_index": "/mnt/analysis/marsoc/pipestances/HC7WVDMXX/BCL_PROCESSOR_PD/HC7WVDMXX/1015.12.13-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/MERGE_FASTQS_FROM_TILES/fork0/join-u37c6ba9993/files/demultiplexed_fastq_path/read-I1_si-CGGAGCAC_lane-002-chunk-001.fastq.gz",
            "subsample_rate": 0.37474273634185645
        },
        {
            "barcode": "/mnt/analysis/marsoc/pipestances/HC7WVDMXX/BCL_PROCESSOR_PD/HC7WVDMXX/1015.12.13-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/MERGE_FASTQS_FROM_TILES/fork0/join-u37c6ba9993/files/demultiplexed_fastq_path/read-I2_si-GACCTATT_lane-001-chunk-001.fastq.gz",
            "barcode_reverse_complement": false,
            "gem_group": 1,
            "read1": "/mnt/analysis/marsoc/pipestances/HC7WVDMXX/BCL_PROCESSOR_PD/HC7WVDMXX/1015.12.13-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/MERGE_FASTQS_FROM_TILES/fork0/join-u37c6ba9993/files/demultiplexed_fastq_path/read-RA_si-GACCTATT_lane-001-chunk-001.fastq.gz",
            "read2": null,
            "read_group": "66333:66333:1:HC7WVDMXX:1",
            "reads_interleaved": true,
            "sample_index": "/mnt/analysis/marsoc/pipestances/HC7WVDMXX/BCL_PROCESSOR_PD/HC7WVDMXX/1015.12.13-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/MERGE_FASTQS_FROM_TILES/fork0/join-u37c6ba9993/files/demultiplexed_fastq_path/read-I1_si-GACCTATT_lane-001-chunk-001.fastq.gz",
            "subsample_rate": 0.37474273634185645
        },
        {
            "barcode": "/mnt/analysis/marsoc/pipestances/HC7WVDMXX/BCL_PROCESSOR_PD/HC7WVDMXX/1015.12.13-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/MERGE_FASTQS_FROM_TILES/fork0/join-u37c6ba9993/files/demultiplexed_fastq_path/read-I2_si-GACCTATT_lane-002-chunk-001.fastq.gz",
            "barcode_reverse_complement": false,
            "gem_group": 1,
            "read1": "/mnt/analysis/marsoc/pipestances/HC7WVDMXX/BCL_PROCESSOR_PD/HC7WVDMXX/1015.12.13-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/MERGE_FASTQS_FROM_TILES/fork0/join-u37c6ba9993/files/demultiplexed_fastq_path/read-RA_si-GACCTATT_lane-002-chunk-001.fastq.gz",
            "read2": null,
            "read_group": "66333:66333:1:HC7WVDMXX:2",
            "reads_interleaved": true,
            "sample_index": "/mnt/analysis/marsoc/pipestances/HC7WVDMXX/BCL_PROCESSOR_PD/HC7WVDMXX/1015.12.13-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/MERGE_FASTQS_FROM_TILES/fork0/join-u37c6ba9993/files/demultiplexed_fastq_path/read-I1_si-GACCTATT_lane-002-chunk-001.fastq.gz",
            "subsample_rate": 0.37474273634185645
        },
        {
            "barcode": "/mnt/analysis/marsoc/pipestances/HC7WVDMXX/BCL_PROCESSOR_PD/HC7WVDMXX/1015.12.13-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/MERGE_FASTQS_FROM_TILES/fork0/join-u37c6ba9993/files/demultiplexed_fastq_path/read-I2_si-ACTTAGGA_lane-001-chunk-001.fastq.gz",
            "barcode_reverse_complement": false,
            "gem_group": 1,
            "read1": "/mnt/analysis/marsoc/pipestances/HC7WVDMXX/BCL_PROCESSOR_PD/HC7WVDMXX/1015.12.13-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/MERGE_FASTQS_FROM_TILES/fork0/join-u37c6ba9993/files/demultiplexed_fastq_path/read-RA_si-ACTTAGGA_lane-001-chunk-001.fastq.gz",
            "read2": null,
            "read_group": "66333:66333:1:HC7WVDMXX:1",
            "reads_interleaved": true,
            "sample_index": "/mnt/analysis/marsoc/pipestances/HC7WVDMXX/BCL_PROCESSOR_PD/HC7WVDMXX/1015.12.13-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/MERGE_FASTQS_FROM_TILES/fork0/join-u37c6ba9993/files/demultiplexed_fastq_path/read-I1_si-ACTTAGGA_lane-001-chunk-001.fastq.gz",
            "subsample_rate": 0.37474273634185645
        },
        {
            "barcode": "/mnt/analysis/marsoc/pipestances/HC7WVDMXX/BCL_PROCESSOR_PD/HC7WVDMXX/1015.12.13-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/MERGE_FASTQS_FROM_TILES/fork0/join-u37c6ba9993/files/demultiplexed_fastq_path/read-I2_si-ACTTAGGA_lane-002-chunk-001.fastq.gz",
            "barcode_reverse_complement": false,
            "gem_group": 1,
            "read1": "/mnt/analysis/marsoc/pipestances/HC7WVDMXX/BCL_PROCESSOR_PD/HC7WVDMXX/1015.12.13-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/MERGE_FASTQS_FROM_TILES/fork0/join-u37c6ba9993/files/demultiplexed_fastq_path/read-RA_si-ACTTAGGA_lane-002-chunk-001.fastq.gz",
            "read2": null,
            "read_group": "66333:66333:1:HC7WVDMXX:2",
            "reads_interleaved": true,
            "sample_index": "/mnt/analysis/marsoc/pipestances/HC7WVDMXX/BCL_PROCESSOR_PD/HC7WVDMXX/1015.12.13-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/MERGE_FASTQS_FROM_TILES/fork0/join-u37c6ba9993/files/demultiplexed_fastq_path/read-I1_si-ACTTAGGA_lane-002-chunk-001.fastq.gz",
            "subsample_rate": 0.37474273634185645
        },
        {
            "barcode": "/mnt/analysis/marsoc/pipestances/HC7WVDMXX/BCL_PROCESSOR_PD/HC7WVDMXX/1015.12.13-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/MERGE_FASTQS_FROM_TILES/fork0/join-u37c6ba9993/files/demultiplexed_fastq_path/read-I2_si-TTAGCTCG_lane-001-chunk-001.fastq.gz",
            "barcode_reverse_complement": false,
            "gem_group": 1,
            "read1": "/mnt/analysis/marsoc/pipestances/HC7WVDMXX/BCL_PROCESSOR_PD/HC7WVDMXX/1015.12.13-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/MERGE_FASTQS_FROM_TILES/fork0/join-u37c6ba9993/files/demultiplexed_fastq_path/read-RA_si-TTAGCTCG_lane-001-chunk-001.fastq.gz",
            "read2": null,
            "read_group": "66333:66333:1:HC7WVDMXX:1",
            "reads_interleaved": true,
            "sample_index": "/mnt/analysis/marsoc/pipestances/HC7WVDMXX/BCL_PROCESSOR_PD/HC7WVDMXX/1015.12.13-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/MERGE_FASTQS_FROM_TILES/fork0/join-u37c6ba9993/files/demultiplexed_fastq_path/read-I1_si-TTAGCTCG_lane-001-chunk-001.fastq.gz",
            "subsample_rate": 0.37474273634185645
        },
        {
            "barcode": "/mnt/analysis/marsoc/pipestances/HC7WVDMXX/BCL_PROCESSOR_PD/HC7WVDMXX/1015.12.13-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/MERGE_FASTQS_FROM_TILES/fork0/join-u37c6ba9993/files/demultiplexed_fastq_path/read-I2_si-TTAGCTCG_lane-002-chunk-001.fastq.gz",
            "barcode_reverse_complement": false,
            "gem_group": 1,
            "read1": "/mnt/analysis/marsoc/pipestances/HC7WVDMXX/BCL_PROCESSOR_PD/HC7WVDMXX/1015.12.13-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/MERGE_FASTQS_FROM_TILES/fork0/join-u37c6ba9993/files/demultiplexed_fastq_path/read-RA_si-TTAGCTCG_lane-002-chunk-001.fastq.gz",
            "read2": null,
            "read_group": "66333:66333:1:HC7WVDMXX:2",
            "reads_interleaved": true,
            "sample_index": "/mnt/analysis/marsoc/pipestances/HC7WVDMXX/BCL_PROCESSOR_PD/HC7WVDMXX/1015.12.13-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/MERGE_FASTQS_FROM_TILES/fork0/join-u37c6ba9993/files/demultiplexed_fastq_path/read-I1_si-TTAGCTCG_lane-002-chunk-001.fastq.gz",
            "subsample_rate": 0.37474273634185645
        }
    ]
    "#;

    const SCDNA_CFG: &str = r#"
    [
        {
            "barcode": null,
            "barcode_reverse_complement": false,
            "bc_in_read": 1,
            "bc_length": 16,
            "gem_group": 1,
            "read1": "/mnt/analysis/marsoc/pipestances/H3TNHDSXX/BCL_PROCESSOR_PD/H3TNHDSXX/1015.12.13-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/MERGE_FASTQS_FROM_TILES/fork0/join-u2388c679b0/files/demultiplexed_fastq_path/read-RA_si-GCGGATAG_lane-002-chunk-001.fastq.gz",
            "read2": null,
            "read_group": "68156:68156:1:unknown_fc:0",
            "reads_interleaved": true,
            "sample_index": "/mnt/analysis/marsoc/pipestances/H3TNHDSXX/BCL_PROCESSOR_PD/H3TNHDSXX/1015.12.13-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/MERGE_FASTQS_FROM_TILES/fork0/join-u2388c679b0/files/demultiplexed_fastq_path/read-I1_si-GCGGATAG_lane-002-chunk-001.fastq.gz",
            "subsample_rate": 1.0
        },
        {
            "barcode": null,
            "barcode_reverse_complement": false,
            "bc_in_read": 1,
            "bc_length": 16,
            "gem_group": 1,
            "read1": "/mnt/analysis/marsoc/pipestances/H3TNHDSXX/BCL_PROCESSOR_PD/H3TNHDSXX/1015.12.13-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/MERGE_FASTQS_FROM_TILES/fork0/join-u2388c679b0/files/demultiplexed_fastq_path/read-RA_si-CATCTCCA_lane-002-chunk-001.fastq.gz",
            "read2": null,
            "read_group": "68156:68156:1:unknown_fc:0",
            "reads_interleaved": true,
            "sample_index": "/mnt/analysis/marsoc/pipestances/H3TNHDSXX/BCL_PROCESSOR_PD/H3TNHDSXX/1015.12.13-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/MERGE_FASTQS_FROM_TILES/fork0/join-u2388c679b0/files/demultiplexed_fastq_path/read-I1_si-CATCTCCA_lane-002-chunk-001.fastq.gz",
            "subsample_rate": 1.0
        },
        {
            "barcode": null,
            "barcode_reverse_complement": false,
            "bc_in_read": 1,
            "bc_length": 16,
            "gem_group": 1,
            "read1": "/mnt/analysis/marsoc/pipestances/H3TNHDSXX/BCL_PROCESSOR_PD/H3TNHDSXX/1015.12.13-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/MERGE_FASTQS_FROM_TILES/fork0/join-u2388c679b0/files/demultiplexed_fastq_path/read-RA_si-AGCAGAGC_lane-002-chunk-001.fastq.gz",
            "read2": null,
            "read_group": "68156:68156:1:unknown_fc:0",
            "reads_interleaved": true,
            "sample_index": "/mnt/analysis/marsoc/pipestances/H3TNHDSXX/BCL_PROCESSOR_PD/H3TNHDSXX/1015.12.13-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/MERGE_FASTQS_FROM_TILES/fork0/join-u2388c679b0/files/demultiplexed_fastq_path/read-I1_si-AGCAGAGC_lane-002-chunk-001.fastq.gz",
            "subsample_rate": 1.0
        },
        {
            "barcode": null,
            "barcode_reverse_complement": false,
            "bc_in_read": 1,
            "bc_length": 16,
            "gem_group": 1,
            "read1": "/mnt/analysis/marsoc/pipestances/H3TNHDSXX/BCL_PROCESSOR_PD/H3TNHDSXX/1015.12.13-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/MERGE_FASTQS_FROM_TILES/fork0/join-u2388c679b0/files/demultiplexed_fastq_path/read-RA_si-TTATCGTT_lane-002-chunk-001.fastq.gz",
            "read2": null,
            "read_group": "68156:68156:1:unknown_fc:0",
            "reads_interleaved": true,
            "sample_index": "/mnt/analysis/marsoc/pipestances/H3TNHDSXX/BCL_PROCESSOR_PD/H3TNHDSXX/1015.12.13-0/BCL_PROCESSOR_PD/BCL_PROCESSOR/MERGE_FASTQS_FROM_TILES/fork0/join-u2388c679b0/files/demultiplexed_fastq_path/read-I1_si-TTATCGTT_lane-002-chunk-001.fastq.gz",
            "subsample_rate": 1.0
        }
    ]
    "#;
}