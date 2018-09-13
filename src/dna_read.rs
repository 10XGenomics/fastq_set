// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

//! ReadPair wrapper object for DNA reads from Linked-Read Genome,  Single-Cell DNA,
//! and Single-Cell ATAC libraries. Provides access to the barcode and allows for dynamic
//! trimming.

use read_pair::{ReadPair, ReadPart, RpRange, WhichRead};
use std::ops::Range;
use {AlignableReadPair, Barcode, FastqProcessor, HasBamTags, HasBarcode, InputFastqs};

/// Specification of a single set of FASTQs and how to interpret the read components.
/// This structure is produced by the SETUP_CHUNKS stage of DNA pipelines
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
}

/// Process raw FASTQ data into DnaRead objects, based on the DnaChunk parameters.
/// This is achieved by implementing `FastqProcessor<DnaRead>`.
#[derive(Clone)]
pub struct DnaProcessor {
    chunk: DnaChunk,
    chunk_id: u16,
}

impl FastqProcessor for DnaProcessor {
    type ReadType = DnaRead;
    /// Convert raw FASTQ data into DnaRead
    fn process_read(&self, read: ReadPair) -> Option<DnaRead> {
        assert!(read.get(WhichRead::R1, ReadPart::Seq).is_some());
        assert!(read.get(WhichRead::R2, ReadPart::Seq).is_some());

        let chunk = &self.chunk;
        // Setup initial (uncorrected) bacode
        let bc_length = chunk.bc_length.unwrap_or(16);

        let bc_range = match chunk.bc_in_read {
            Some(1) => RpRange::new(WhichRead::R1, 0, Some(bc_length)),
            None => RpRange::new(WhichRead::I2, 0, chunk.bc_length),
            _ => panic!("unsupported barcode location"),
        };

        // Snip out barcode
        let barcode = Barcode::new(
            chunk.gem_group,
            read.get_range(&bc_range, ReadPart::Seq).unwrap(),
            true,
        );

        // FIXME -- do cutadapt trimming here?

        Some(DnaRead {
            data: read,
            barcode,
            bc_range,
            chunk_id: self.chunk_id,

            // FIXME -- get the right trim length
            trim_r1: 7,
            trim_r2: 0,
        })
    }

    fn fastq_files(&self) -> InputFastqs {
        InputFastqs {
            r1: self.chunk.read1.clone(),
            r2: self.chunk.read2.clone(),
            i1: self.chunk.sample_index.clone(),
            i2: self.chunk.barcode.clone(),
            r1_interleaved: self.chunk.reads_interleaved,
        }
    }

    fn bc_subsample_rate(&self) -> f64 {
        1.0
    }

    fn read_subsample_rate(&self) -> f64 {
        self.chunk.subsample_rate
    }

    fn gem_group(&self) -> u16 {
        self.chunk.gem_group
    }
}

/// Represents a GEM-barcoded DNA read, with a barcode at the start of R1 or in an index read,
/// and possibly some bases trimmed the the start of R1 and R2.
#[derive(Serialize, Deserialize, Eq, PartialEq)]
pub struct DnaRead {
    data: ReadPair,
    barcode: Barcode,
    bc_range: RpRange,
    trim_r1: u8,
    trim_r2: u8,
    chunk_id: u16,
}

impl HasBarcode for DnaRead {
    fn barcode(&self) -> &Barcode {
        &self.barcode
    }

    fn set_barcode(&mut self, barcode: Barcode) {
        self.barcode = barcode
    }

    fn barcode_qual(&self) -> &[u8] {
        self.raw_bc_seq()
    }
}

impl HasBamTags for DnaRead {
    fn tags(&self) -> Vec<([u8; 2], &[u8])> {
        vec![
            (*b"RX", self.raw_bc_seq()),
            (*b"QX", self.raw_bc_qual()),
            (*b"TR", self.r1_trim_seq()),
            (*b"TQ", self.r1_trim_qual()),
        ]
    }
}

impl DnaRead {
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
        self.data.get_range(&self.bc_range, ReadPart::Seq).unwrap()
    }

    /// Raw barcode QVs
    pub fn raw_bc_qual(&self) -> &[u8] {
        self.data.get_range(&self.bc_range, ReadPart::Qual).unwrap()
    }

    #[inline]
    pub fn r1_trim_range(&self) -> Range<usize> {
        if self.bc_range.read() == WhichRead::R1 {
            let bcr = self.bc_range;
            let start = bcr.offset() + bcr.len().unwrap_or(0);

            start..start + self.trim_r1 as usize
        } else {
            0..self.trim_r1 as usize
        }
    }

    /// Bases trimmed after the 10x BC, before the start of bases used from R1
    pub fn r1_trim_seq(&self) -> &[u8] {
        let rng = self.r1_trim_range();
        &self.r1_seq_raw()[rng]
    }

    /// QVs trimmed after the 10x BC, before the start of bases used from R1
    pub fn r1_trim_qual(&self) -> &[u8] {
        let rng = self.r1_trim_range();
        &self.r1_qual_raw()[rng]
    }

    /// Usable R1 bases after removal of BC and trimming
    pub fn r1_seq(&self) -> &[u8] {
        let rng = self.r1_trim_range();
        &self.r1_seq_raw()[rng.end..]
    }

    /// Usable R1 bases after removal of BC and trimming
    pub fn r1_qual(&self) -> &[u8] {
        let rng = self.r1_trim_range();
        &self.r1_qual_raw()[rng.end..]
    }
}

impl AlignableReadPair for DnaRead {
    fn header(&self) -> &[u8] {
        self.header()
    }

    fn alignable_sequence(&self) -> (&[u8], &[u8]) {
        (self.r1_seq(), self.r2_seq())
    }

    fn alignable_quals(&self) -> (&[u8], &[u8]) {
        (self.r1_qual(), self.r2_qual())
    }
}

#[cfg(test)]
mod test_dna_cfg {
    use super::*;
    use serde_json;

    fn load_dna_chunk_def(chunk_json: &str) -> Vec<DnaChunk> {
        serde_json::from_str(chunk_json).unwrap()
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
        let chunks = load_dna_chunk_def(ATAC_CFG_TEST);
        println!("{:?}", chunks);

        let mut procs = Vec::new();
        for (idx, chunk) in chunks.into_iter().enumerate() {
            let prc = DnaProcessor {
                chunk: chunk,
                chunk_id: idx as u16,
            };

            procs.push(prc);
        }

        // let bc_sort_results =
            // ::bc_sort::barcode_sort_workflow(procs, "test", "test/10K-agora-dev.txt").unwrap();
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
            "subsample_rate": 1.0
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
            "subsample_rate": 1.0
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
