// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

//! ReadPair wrapper object for RNA reads from Single Cell 3' nad Single Cell 5' / VDJ ibraries.
//! Provides access to the barcode and allows for dynamic trimming.

use adapter_trimmer::{intersect_ranges, ReadAdapterCatalog, AdapterTrimmer};
use read_pair::{ReadPair, ReadPart, RpRange, WhichRead};
use fxhash::FxHashMap;
use std::ops::Range;
use WhichEnd;
use {Barcode, FastqProcessor, HasBarcode, HasSampleIndex, InputFastqs, Umi};

#[derive(Serialize, Deserialize, Clone, PartialEq, Eq, Debug)]
/// Define a chemistry supported by our RNA products.
///
/// A chemistry tells you where & how to look for various read components
/// (cell barcode, cDNA sequence, UMI, sample index) from the FASTQ
/// cluster data. It is convenient to specify the chemistry as a JSON
/// file and let `ChemistryDef` deserialize it.
///
/// As an example, consider the Single Cell V(D)J read layout below:
/// ![Plot](../../../../doc-media/fastq/scvdj_chem.png)
/// The corresponding chemistry definition would be:
///
/// - Barcode is present in the first 16 bases of Read 1.
/// This translates to a `read_type` "R1", `read_offset` 0
/// and `read_length` 16. Valid options for `read_type` are
/// "R1", "R2", "I1", "I2"
/// ``` text
/// {
///     "barcode_read_length": 16,
///     "barcode_read_offset": 0,
///     "barcode_read_type": "R1",
///     "barcode_whitelist": "737K-august-2016",
/// ```
/// - Description and name for the chemistry
/// ``` text
///     "description": "Single Cell V(D)J",
///     "name": "SCVDJ",
/// ```
/// - Specify the `endedness` of the product. This would be `three_prime` for Single Cell
/// 3' Gene expression and `five_prime` for VDJ and 5' Gene expression
/// ``` text
///     "endedness": "five_prime",
/// ```
/// - Filename conventions used by `bcl2fastq` and `bcl_processor`
/// ``` text
///     "read_type_to_bcl2fastq_filename": {
///         "I1": "I1",
///         "I2": null,
///         "R1": "R1",
///         "R2": "R2"
///     },
///     "read_type_to_bcl_processor_filename": {
///         "I1": "I1",
///         "I2": null,
///         "R1": "RA",
///         "R2": null
///     },
/// ```
/// - Every RNA product will have cDNA sequences in
/// one or both of read1 and read2. Hence an `rna_read`
/// definition is necessary for each chemistry. `rna_read2`
/// is optional. For V(D)J, `rna_read` would be the bases in read1
/// beyond the spacer sequence and `rna_read2` would be all the
/// bases in read2. It is not necessary that `rna_read` correspond to
/// read1. For example, in Single Cell 5' R2-only mode, `rna_read`
/// would be all of read2 and `rna_read2` would be empty.
/// ``` text
///     "rna_read_length": null,
///     "rna_read_offset": 41,
///     "rna_read_type": "R1",
/// ```
/// - Optionally specify `rna_read2`
/// ``` text
///     "rna_read2_length": null,
///     "rna_read2_offset": 0,
///     "rna_read2_type": "R2",
/// ```
/// - Sample Index read definition
/// ``` text
///     "si_read_length": null,
///     "si_read_offset": 0,
///     "si_read_type": "I1",
/// ```
/// - Strandedness is `+` when the rna_read and the transcript are
/// expected to be in the same orientation and `-` otherwise.
/// ``` text
///     "strandedness": "+",
/// ```
/// - UMI is present in the bases 16-25 of read1
/// ``` text
///     "umi_read_length": 10,
///     "umi_read_offset": 16,
///     "umi_read_type": "R1"
/// }
/// ```
pub struct ChemistryDef {
    barcode_read_length: usize,
    barcode_read_offset: usize,
    barcode_read_type: WhichRead,
    barcode_whitelist: String,
    description: String,
    endedness: Option<WhichEnd>,
    name: String,
    read_type_to_bcl2fastq_filename: FxHashMap<WhichRead, Option<String>>,
    read_type_to_bcl_processor_filename: FxHashMap<WhichRead, Option<String>>,
    rna_read_length: Option<usize>,
    rna_read_offset: usize,
    rna_read_type: WhichRead,
    rna_read2_length: Option<usize>,
    rna_read2_offset: Option<usize>,
    rna_read2_type: Option<WhichRead>,
    si_read_length: Option<usize>,
    si_read_offset: usize,
    si_read_type: WhichRead,
    strandedness: String,
    umi_read_length: usize,
    umi_read_offset: usize,
    umi_read_type: WhichRead,
}

impl ChemistryDef {
    pub fn is_paired_end(&self) -> bool {
        self.rna_read2_type.is_some()
    }
}

#[derive(Serialize, Deserialize, Clone, PartialEq, Debug)]
struct ReadChunks {
    #[serde(rename = "R1")]
    r1: String,
    #[serde(rename = "R2")]
    r2: Option<String>,
    #[serde(rename = "I1")]
    i1: Option<String>,
    #[serde(rename = "I2")]
    i2: Option<String>,
}

#[derive(Serialize, Deserialize, Clone, PartialEq, Debug)]
pub struct RnaChunk {
    chemistry: ChemistryDef,
    gem_group: u16,
    read_chunks: ReadChunks,
    read_group: String,
    reads_interleaved: bool,
    subsample_rate: Option<f64>,
    #[serde(default)]
    read_lengths: FxHashMap<WhichRead, usize>,
}

impl RnaChunk {
    pub fn subsample_rate(&mut self, value: f64) -> &mut Self {
        self.subsample_rate = Some(value);
        self
    }
    pub fn illumina_r1_trim_length(&mut self, value: usize) -> &mut Self {
        self.read_lengths.insert(WhichRead::R1, value);
        self
    }
    pub fn illumina_r2_trim_length(&mut self, value: usize) -> &mut Self {
        self.read_lengths.insert(WhichRead::R2, value);
        self
    }
}

impl FastqProcessor for RnaChunk {
    type ReadType = RnaRead;
    fn process_read(&self, read: ReadPair) -> Option<RnaRead> {
        let chem = &self.chemistry;

        let bc_range = RpRange::new(
            chem.barcode_read_type,
            chem.barcode_read_offset,
            Some(chem.barcode_read_length),
        );
        let umi_range = RpRange::new(
            chem.umi_read_type,
            chem.umi_read_offset,
            Some(chem.umi_read_length),
        );
        let r1_range = {
            let read_type = chem.rna_read_type;
            let mut range = RpRange::new(
                read_type,
                chem.rna_read_offset,
                chem.rna_read_length,
            );
            if let Some(len) = self.read_lengths.get(&read_type) {
                range.intersect(RpRange::new(read_type, 0, Some(*len)))
            }
            range
        };

        let r2_range = match chem.rna_read2_type {
            Some(read_type) => {
                let mut range = RpRange::new(
                    read_type,
                    chem.rna_read2_offset.unwrap(),
                    chem.rna_read2_length,
                );
                if let Some(len) = self.read_lengths.get(&read_type) {
                    range.intersect(RpRange::new(read_type, 0, Some(*len)))
                }
                Some(range)
            },
            None => None,
        };

        let barcode = Barcode::new(
            self.gem_group,
            read.get_range(&bc_range, ReadPart::Seq).unwrap(),
            true,
        );
        let umi = Umi::new(read.get_range(&umi_range, ReadPart::Seq).unwrap());

        Some(RnaRead {
            read,
            barcode,
            umi,
            bc_range,
            umi_range,
            r1_range,
            r2_range,
        })
    }

    fn fastq_files(&self) -> InputFastqs {
        // Make sure that either 
        // - r2 is None and interleaved
        // - r2 is Some and not interleaved
        match self.read_chunks.r2 {
            Some(_) => assert!(!self.reads_interleaved),
            None => assert!(self.reads_interleaved),
        }
        
        InputFastqs {
            r1: self.read_chunks.r1.clone(),
            r2: self.read_chunks.r2.clone(),
            i1: self.read_chunks.i1.clone(),
            i2: self.read_chunks.i2.clone(),
            r1_interleaved: self.reads_interleaved,
        }
    }
    fn bc_subsample_rate(&self) -> f64 {
        1.0
    }
    fn read_subsample_rate(&self) -> f64 {
        self.subsample_rate.unwrap_or(1.0)
    }

    fn gem_group(&self) -> u16 {
        self.gem_group
    }
}

#[derive(Serialize, Deserialize, Debug, Eq, PartialEq, Clone)]
pub struct RnaRead {
    read: ReadPair,
    barcode: Barcode,
    umi: Umi,
    bc_range: RpRange,
    umi_range: RpRange,
    r1_range: RpRange,
    r2_range: Option<RpRange>,
}

impl HasBarcode for RnaRead {
    fn barcode(&self) -> &Barcode {
        &self.barcode
    }

    fn set_barcode(&mut self, barcode: Barcode) {
        self.barcode = barcode;
    }

    fn barcode_qual(&self) -> &[u8] {
        self.raw_bc_qual()
    }

    fn raw_bc_seq(&self) -> &[u8] {
        self.read.get_range(&self.bc_range, ReadPart::Seq).unwrap()
    }

    fn raw_bc_qual(&self) -> &[u8] {
        self.read.get_range(&self.bc_range, ReadPart::Qual).unwrap()
    }
}

impl HasSampleIndex for RnaRead {
    fn si_seq(&self) -> Option<&[u8]> {
        self.read.get(WhichRead::I1, ReadPart::Seq)
    }

    fn si_qual(&self) -> Option<&[u8]> {
        self.read.get(WhichRead::I1, ReadPart::Qual)
    }
}

impl RnaRead {
    pub fn bc_range(&self) -> &RpRange {
        &self.bc_range
    }
    pub fn umi_range(&self) -> &RpRange {
        &self.umi_range
    }
    pub fn readpair(&self) -> &ReadPair {
        &self.read
    }

    pub fn is_paired_end(&self) -> bool {
        self.r2_range.is_some()
    }

    /// FASTQ read header
    pub fn header(&self) -> &[u8] {
        self.read.get(WhichRead::R1, ReadPart::Header).unwrap()
    }

    /// Full raw R1 sequence
    pub fn raw_illumina_read1_seq(&self) -> &[u8] {
        self.read.get(WhichRead::R1, ReadPart::Seq).unwrap()
    }

    /// Full raw R1 QVs
    pub fn raw_illumina_read1_qual(&self) -> &[u8] {
        self.read.get(WhichRead::R1, ReadPart::Qual).unwrap()
    }

    /// Full R2 sequence
    pub fn raw_illumina_read2_seq(&self) -> &[u8] {
        self.read.get(WhichRead::R2, ReadPart::Seq).unwrap()
    }

    /// Full R2 QVs
    pub fn raw_illumina_read2_qual(&self) -> &[u8] {
        self.read.get(WhichRead::R2, ReadPart::Qual).unwrap()
    }

    /// Raw, uncorrected barcode sequence
    pub fn raw_umi_seq(&self) -> &[u8] {
        self.read.get_range(&self.umi_range, ReadPart::Seq).unwrap()
    }

    /// Raw barcode QVs
    pub fn raw_umi_qual(&self) -> &[u8] {
        self.read
            .get_range(&self.umi_range, ReadPart::Qual)
            .unwrap()
    }

    /// Usable R1 bases after removal of BC and trimming
    pub fn r1_seq(&self) -> &[u8] {
        self.read.get_range(&self.r1_range, ReadPart::Seq).unwrap()
    }

    /// Usable R1 bases after removal of BC and trimming
    pub fn r1_qual(&self) -> &[u8] {
        self.read.get_range(&self.r1_range, ReadPart::Qual).unwrap()
    }

    /// Usable R2 bases after removal of BC and trimming
    pub fn r2_seq(&self) -> Option<&[u8]> {
        if let Some(range) = self.r2_range {
            self.read.get_range(&range, ReadPart::Seq)
        } else {
            None
        }
    }

    /// Usable R2 bases after removal of BC and trimming
    pub fn r2_qual(&self) -> Option<&[u8]> {
        if let Some(range) = self.r2_range {
            self.read.get_range(&range, ReadPart::Qual)
        } else {
            None
        }
    }

    fn trim_adapters_helper<'a>(
        &mut self,
        read_range: RpRange,
        trimmers: &mut Vec<AdapterTrimmer<'a>>,
    ) -> (Range<usize>, FxHashMap<String, RpRange>) {
        let seq = self.read.get_range(&read_range, ReadPart::Seq).unwrap();
        let trim_results: Vec<_> = trimmers
            .iter_mut()
            .filter_map(|t| {
                t.find(seq)
                    .map(|trim_result| (t.adapter.name.clone(), trim_result))
            })
            .collect();
        let shrink_range = trim_results.iter().fold(0..seq.len(), |acc, (_, x)| {
            intersect_ranges(&acc, &x.retain_range)
        });
        let adapter_pos = trim_results.into_iter().fold(FxHashMap::default(), |mut acc, (name, trim_result)| {
            let mut this_range = read_range.clone();
            this_range.shrink(&trim_result.adapter_range);
            acc.insert(name, this_range);
            acc
        });
        (shrink_range, adapter_pos)
    }

    pub fn trim_adapters<'a>(
        &mut self,
        adapter_catalog: &'a mut ReadAdapterCatalog<'a>,
    ) -> FxHashMap<String, RpRange> {

        let mut result = FxHashMap::default();
        // Trim r1
        {
            let ad_trimmers = adapter_catalog.get_mut_trimmers(self.r1_range.read());
            let range = self.r1_range; // Creates a copy
            let (shrink_range, adapter_pos) = self.trim_adapters_helper(range, ad_trimmers);
            result.extend(adapter_pos);
            self.r1_range.shrink(&shrink_range);
        }

        // Trim r2
        {
            if let Some(range) = self.r2_range {
                let ad_trimmers = adapter_catalog.get_mut_trimmers(range.read());
                let (shrink_range, adapter_pos) = self.trim_adapters_helper(range, ad_trimmers);
                result.extend(adapter_pos);
                self.r2_range.as_mut().map(|x| x.shrink(&shrink_range));
            }
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json;

    #[test]
    fn test_vdj_cfg() {
        use std::fs::File;
        let c: Vec<RnaChunk> = serde_json::from_reader(File::open("tests/vdj_rna_chunk.json").unwrap()).unwrap();
        println!("{:?}", c);
    }

    #[test]
    fn test_gex_cfg() {
        use std::fs::File;
        let c: Vec<RnaChunk> = serde_json::from_reader(File::open("tests/gex_rna_chunk.json").unwrap()).unwrap();
        println!("{:?}", c);
    }
}
