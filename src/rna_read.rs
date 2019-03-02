// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

//! ReadPair wrapper object for RNA reads from Single Cell 3' nad Single Cell 5' / VDJ ibraries.
//! Provides access to the barcode and allows for dynamic trimming.

use adapter_trimmer::{intersect_ranges, ReadAdapterCatalog, AdapterTrimmer};
use read_pair::{ReadPair, ReadPart, RpRange, WhichRead};
use fxhash::FxHashMap;
use std::ops::Range;
use WhichEnd;
use {Barcode, FastqProcessor, HasBarcode, HasSampleIndex, InputFastqs, Umi};
use std::cmp::min;

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

/// `RnaChunk` represents a chunk of reads from any of our RNA products.
/// This is typically created by the `SETUP_CHUNKS` stage in our pipelines.
/// This struct knows how to interpret the raw fastq data using the `chemistry`.
/// `RnaChunk` is a `FastqProcessor`, meaning it can process the raw reads and
/// create `RnaRead`, which resolves various regions in the read such as BC, UMI etc.
/// You can also specify a subsample rate and set R1/R2 lengths.
/// 
/// # Example
/// See the tests
/// 
/// # Tests
/// * `test_rna_chunk_processor_interleaved_sc_vdj` - Test that the `RnaRead`
/// produced by the `RnaChunk` has the expected barcode, umi, read1, read2 for
/// VDJ chemistry with interleaved reads.
/// * `test_rna_chunk_processor_interleaved_sc_vdj_r2` - Test that the `RnaRead`
/// produced by the `RnaChunk` has the expected barcode, umi, read1, read2 for
/// VDJ R2-only chemistry with interleaved reads.
/// * `test_rna_chunk_processor_non_interleaved_sc_vdj` - Test that the `RnaRead`
/// produced by the `RnaChunk` has the expected barcode, umi, read1, read2 for
/// VDJ chemistry with non-interleaved reads.
/// * `test_rna_chunk_processor_non_interleaved_sc_vdj_r2` - Test that the `RnaRead`
/// produced by the `RnaChunk` has the expected barcode, umi, read1, read2 for
/// VDJ R2-only chemistry with non-interleaved reads.
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

    /// Create a new `RnaChunk`.
    pub fn new(chemistry: ChemistryDef, gem_group: u16, fastqs: InputFastqs, read_group: String) -> RnaChunk {
        RnaChunk {
            chemistry,
            gem_group,
            read_chunks: ReadChunks {
                r1: fastqs.r1,
                r2: fastqs.r2,
                i1: fastqs.i1,
                i2: fastqs.i2,
            },
            read_group,
            reads_interleaved: fastqs.r1_interleaved,
            subsample_rate: None,
            read_lengths: FxHashMap::default(),
        }
    }

    /// Set the subsample rate
    /// 
    /// # Test
    /// * `test_rna_chunk_subsample()` - Make sure we get roughly as many
    /// reads as expected after subsampling
    pub fn subsample_rate(&mut self, value: f64) -> &mut Self {
        self.subsample_rate = Some(value);
        self
    }
    /// Set the length to hard trim the read1 in the input fastq.
    /// 
    /// # Test
    /// * `prop_test_rna_chunk_trim()` - Make sure that trimming works as
    /// expected for arbitrary inputs
    pub fn illumina_r1_trim_length(&mut self, value: usize) -> &mut Self {
        self.read_lengths.insert(WhichRead::R1, value);
        self
    }
    /// Set the length to hard trim the read2 in the input fastq
    /// 
    /// # Test
    /// * `prop_test_rna_chunk_trim()` - Make sure that trimming works as
    /// expected for arbitrary inputs
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
                let trim_len = min(*len, read.len(read_type).unwrap_or(0));
                range.intersect(RpRange::new(read_type, 0, Some(trim_len)))
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
                    let trim_len = min(*len, read.len(read_type).unwrap_or(0));
                    range.intersect(RpRange::new(read_type, 0, Some(trim_len)))
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

    /// Get the fastq files
    /// 
    /// # Panics
    /// * If interleaved and r2 is not None
    /// * If not interleaved and r2 is None
    /// 
    /// # Tests
    /// * `test_rna_chunk_fastq_panic_1()` - Make sure that we panic if 
    /// interleaved and r2 is not None
    /// * `test_rna_chunk_fastq_panic_2()` - Make sure that we panic if 
    /// not interleaved and r2 is None
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

    /// Raw, uncorrected UMI sequence
    pub fn raw_umi_seq(&self) -> &[u8] {
        self.read.get_range(&self.umi_range, ReadPart::Seq).unwrap()
    }

    /// Raw UMI QVs
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

    /// Given an `RpRange` and a list of adapter trimmers, this function 
    /// searches for all the adapters within the `RpRange` and returns a
    /// `Range<usize>`, which you need to shrink the `RpRange` to.
    /// It also returns a hashmap from the adapter name to the `RpRange`
    /// where the adapter was found. This is useful for book-keeping and 
    /// computing metrics.
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

    /// Perform adapter trimming on the `RnaRead` and return the positions
    /// of the adapters found.
    /// 
    /// # Inputs
    /// * `adapter_catalog`: Packages all the adapter trimmers
    /// 
    /// # Output
    /// * `FxHashMap<String, RpRange>` - where the key is the name of the adapter
    /// and values is the `RpRange` where the adapter is found. Clearly, the adapters
    /// which are not present in the read will not be part of the output. The `ReadAdapterCatalog`
    /// guarantees that no two adapters share the same name, so there is no confusion.ReadPart
    /// 
    /// # Test
    /// * `test_rna_read_adapter_trim()`: Test that we can trim reads consistent with cutadapt.
    pub fn trim_adapters<'a>(
        &mut self,
        adapter_catalog: &mut ReadAdapterCatalog<'a>,
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
    use std::fs::File;
    use std::path::Path;
    use proptest::arbitrary::any;

    #[test]
    fn test_vdj_cfg() {
        let c: Vec<RnaChunk> = serde_json::from_reader(File::open("tests/vdj_rna_chunk.json").unwrap()).unwrap();
        println!("{:?}", c);
    }

    #[test]
    fn test_gex_cfg() {
        let c: Vec<RnaChunk> = serde_json::from_reader(File::open("tests/gex_rna_chunk.json").unwrap()).unwrap();
        println!("{:?}", c);
    }

    fn get_vdj_chemistry() -> ChemistryDef {
        serde_json::from_reader(File::open("tests/rna_read/sc_vdj_chemistry.json").unwrap()).unwrap()
    }

    fn get_vdj_r2_chemistry() -> ChemistryDef {
        serde_json::from_reader(File::open("tests/rna_read/sc_vdj_r2_chemistry.json").unwrap()).unwrap()
    }

    #[test]
    #[should_panic]
    fn test_rna_chunk_fastq_panic_1() {
        let chunk = RnaChunk {
            chemistry: get_vdj_chemistry(),
            gem_group: 1,
            read_chunks: ReadChunks{
                r1: "tests/rna_read/r1_1k.fasta.lz4".into(),
                r2: Some("tests/rna_read/r2_1k.fasta.lz4".into()),
                i1: None,
                i2: None,
            },
            read_group: "Blah".into(),
            reads_interleaved: true,
            subsample_rate: None,
            read_lengths: FxHashMap::default(),
        };
        let _ = chunk.fastq_files();
    }

    #[test]
    #[should_panic]
    fn test_rna_chunk_fastq_panic_2() {
        let chunk = RnaChunk {
            chemistry: get_vdj_chemistry(),
            gem_group: 1,
            read_chunks: ReadChunks{
                r1: "tests/rna_read/r1_1k.fasta.lz4".into(),
                r2: None,
                i1: None,
                i2: None,
            },
            read_group: "Blah".into(),
            reads_interleaved: false,
            subsample_rate: None,
            read_lengths: FxHashMap::default(),
        };
        let _ = chunk.fastq_files();
    }

    fn load_expected(fname: impl AsRef<Path>) -> Vec<Vec<u8>> {
        let expected: Vec<String> = serde_json::from_reader(File::open(fname).unwrap()).unwrap();
        expected.into_iter().map(|x| x.as_bytes().to_vec()).collect()
    }

    #[test]
    fn test_rna_chunk_processor_interleaved_sc_vdj() {
        let chunk = RnaChunk {
            chemistry: get_vdj_chemistry(),
            gem_group: 1,
            read_chunks: ReadChunks{
                r1: "tests/rna_read/interleaved_2k.fastq.lz4".into(),
                r2: None,
                i1: None,
                i2: None,
            },
            read_group: "Blah".into(),
            reads_interleaved: true,
            subsample_rate: None,
            read_lengths: FxHashMap::default(),
        };
        let expected_bc = load_expected("tests/rna_read/bc.json");
        let expected_bc_qual = load_expected("tests/rna_read/bc_qual.json");
        let expected_umi = load_expected("tests/rna_read/umi.json");
        let expected_umi_qual = load_expected("tests/rna_read/umi_qual.json");
        let expected_r1 = load_expected("tests/rna_read/rna_r1.json");
        let expected_r1_qual = load_expected("tests/rna_read/rna_r1_qual.json");
        let expected_r2 = load_expected("tests/rna_read/rna_r2.json");
        let expected_r2_qual = load_expected("tests/rna_read/rna_r2_qual.json");

        for (i, rna_read_result) in chunk.iter().unwrap().enumerate() {
            let rna_read = rna_read_result.unwrap();

            // Check that the barcode data is correct
            assert_eq!(rna_read.bc_range, RpRange::new(WhichRead::R1, 0, Some(16)));
            assert_eq!(rna_read.barcode, Barcode::new(1, &expected_bc[i], true));
            assert_eq!(rna_read.raw_bc_seq(), expected_bc[i].as_slice());
            assert_eq!(rna_read.raw_bc_qual(), expected_bc_qual[i].as_slice());

            // Check that the UMI data is correct
            assert_eq!(rna_read.umi_range, RpRange::new(WhichRead::R1, 16, Some(10)));
            assert_eq!(rna_read.umi, Umi::new(&expected_umi[i]));
            assert_eq!(rna_read.raw_umi_seq(), expected_umi[i].as_slice());
            assert_eq!(rna_read.raw_umi_qual(), expected_umi_qual[i].as_slice());

            // Check that the R1 data is correct
            assert_eq!(rna_read.r1_range, RpRange::new(WhichRead::R1, 41, None));
            assert_eq!(rna_read.r1_seq(), expected_r1[i].as_slice());
            assert_eq!(rna_read.r1_qual(), expected_r1_qual[i].as_slice());

            // Check that the R2 data is correct
            assert_eq!(rna_read.r2_range, Some(RpRange::new(WhichRead::R2, 0, None)));
            assert_eq!(rna_read.r2_seq().unwrap(), expected_r2[i].as_slice());
            assert_eq!(rna_read.r2_qual().unwrap(), expected_r2_qual[i].as_slice());
        }
    }

    #[test]
    fn test_rna_chunk_processor_interleaved_sc_vdj_r2() {
        let chunk = RnaChunk {
            chemistry: get_vdj_r2_chemistry(),
            gem_group: 1,
            read_chunks: ReadChunks{
                r1: "tests/rna_read/interleaved_2k.fastq.lz4".into(),
                r2: None,
                i1: None,
                i2: None,
            },
            read_group: "Blah".into(),
            reads_interleaved: true,
            subsample_rate: None,
            read_lengths: FxHashMap::default(),
        };
        let expected_bc = load_expected("tests/rna_read/bc.json");
        let expected_bc_qual = load_expected("tests/rna_read/bc_qual.json");
        let expected_umi = load_expected("tests/rna_read/umi.json");
        let expected_umi_qual = load_expected("tests/rna_read/umi_qual.json");
        let expected_r1 = load_expected("tests/rna_read/rna_r2.json");
        let expected_r1_qual = load_expected("tests/rna_read/rna_r2_qual.json");

        for (i, rna_read_result) in chunk.iter().unwrap().enumerate() {
            let rna_read = rna_read_result.unwrap();

            // Check that the barcode data is correct
            assert_eq!(rna_read.bc_range, RpRange::new(WhichRead::R1, 0, Some(16)));
            assert_eq!(rna_read.barcode, Barcode::new(1, &expected_bc[i], true));
            assert_eq!(rna_read.raw_bc_seq(), expected_bc[i].as_slice());
            assert_eq!(rna_read.raw_bc_qual(), expected_bc_qual[i].as_slice());

            // Check that the UMI data is correct
            assert_eq!(rna_read.umi_range, RpRange::new(WhichRead::R1, 16, Some(10)));
            assert_eq!(rna_read.umi, Umi::new(&expected_umi[i]));
            assert_eq!(rna_read.raw_umi_seq(), expected_umi[i].as_slice());
            assert_eq!(rna_read.raw_umi_qual(), expected_umi_qual[i].as_slice());

            // Check that the R1 data is correct
            assert_eq!(rna_read.r1_range, RpRange::new(WhichRead::R2, 0, None));
            assert_eq!(rna_read.r1_seq(), expected_r1[i].as_slice());
            assert_eq!(rna_read.r1_qual(), expected_r1_qual[i].as_slice());

            // Check that the R2 data is correct
            assert_eq!(rna_read.r2_range, None);
            assert_eq!(rna_read.r2_seq(), None);
            assert_eq!(rna_read.r2_qual(), None);
        }
    }

    #[test]
    fn test_rna_chunk_processor_non_interleaved_sc_vdj() {
        let chunk = RnaChunk {
            chemistry: get_vdj_chemistry(),
            gem_group: 1,
            read_chunks: ReadChunks{
                r1: "tests/rna_read/r1_2k.fastq.lz4".into(),
                r2: Some("tests/rna_read/r2_2k.fastq.lz4".into()),
                i1: None,
                i2: None,
            },
            read_group: "Blah".into(),
            reads_interleaved: false,
            subsample_rate: None,
            read_lengths: FxHashMap::default(),
        };
        let expected_bc = load_expected("tests/rna_read/bc.json");
        let expected_bc_qual = load_expected("tests/rna_read/bc_qual.json");
        let expected_umi = load_expected("tests/rna_read/umi.json");
        let expected_umi_qual = load_expected("tests/rna_read/umi_qual.json");
        let expected_r1 = load_expected("tests/rna_read/rna_r1.json");
        let expected_r1_qual = load_expected("tests/rna_read/rna_r1_qual.json");
        let expected_r2 = load_expected("tests/rna_read/rna_r2.json");
        let expected_r2_qual = load_expected("tests/rna_read/rna_r2_qual.json");

        for (i, rna_read_result) in chunk.iter().unwrap().enumerate() {
            let rna_read = rna_read_result.unwrap();

            // Check that the barcode data is correct
            assert_eq!(rna_read.bc_range, RpRange::new(WhichRead::R1, 0, Some(16)));
            assert_eq!(rna_read.barcode, Barcode::new(1, &expected_bc[i], true));
            assert_eq!(rna_read.raw_bc_seq(), expected_bc[i].as_slice());
            assert_eq!(rna_read.raw_bc_qual(), expected_bc_qual[i].as_slice());

            // Check that the UMI data is correct
            assert_eq!(rna_read.umi_range, RpRange::new(WhichRead::R1, 16, Some(10)));
            assert_eq!(rna_read.umi, Umi::new(&expected_umi[i]));
            assert_eq!(rna_read.raw_umi_seq(), expected_umi[i].as_slice());
            assert_eq!(rna_read.raw_umi_qual(), expected_umi_qual[i].as_slice());

            // Check that the R1 data is correct
            assert_eq!(rna_read.r1_range, RpRange::new(WhichRead::R1, 41, None));
            assert_eq!(rna_read.r1_seq(), expected_r1[i].as_slice());
            assert_eq!(rna_read.r1_qual(), expected_r1_qual[i].as_slice());

            // Check that the R2 data is correct
            assert_eq!(rna_read.r2_range, Some(RpRange::new(WhichRead::R2, 0, None)));
            assert_eq!(rna_read.r2_seq().unwrap(), expected_r2[i].as_slice());
            assert_eq!(rna_read.r2_qual().unwrap(), expected_r2_qual[i].as_slice());
        }
    }

    #[test]
    fn test_rna_chunk_processor_non_interleaved_sc_vdj_r2() {
        let chunk = RnaChunk {
            chemistry: get_vdj_r2_chemistry(),
            gem_group: 1,
            read_chunks: ReadChunks{
                r1: "tests/rna_read/r1_2k.fastq.lz4".into(),
                r2: Some("tests/rna_read/r2_2k.fastq.lz4".into()),
                i1: None,
                i2: None,
            },
            read_group: "Blah".into(),
            reads_interleaved: false,
            subsample_rate: None,
            read_lengths: FxHashMap::default(),
        };
        let expected_bc = load_expected("tests/rna_read/bc.json");
        let expected_bc_qual = load_expected("tests/rna_read/bc_qual.json");
        let expected_umi = load_expected("tests/rna_read/umi.json");
        let expected_umi_qual = load_expected("tests/rna_read/umi_qual.json");
        let expected_r1 = load_expected("tests/rna_read/rna_r2.json");
        let expected_r1_qual = load_expected("tests/rna_read/rna_r2_qual.json");

        for (i, rna_read_result) in chunk.iter().unwrap().enumerate() {
            let rna_read = rna_read_result.unwrap();

            // Check that the barcode data is correct
            assert_eq!(rna_read.bc_range, RpRange::new(WhichRead::R1, 0, Some(16)));
            assert_eq!(rna_read.barcode, Barcode::new(1, &expected_bc[i], true));
            assert_eq!(rna_read.raw_bc_seq(), expected_bc[i].as_slice());
            assert_eq!(rna_read.raw_bc_qual(), expected_bc_qual[i].as_slice());

            // Check that the UMI data is correct
            assert_eq!(rna_read.umi_range, RpRange::new(WhichRead::R1, 16, Some(10)));
            assert_eq!(rna_read.umi, Umi::new(&expected_umi[i]));
            assert_eq!(rna_read.raw_umi_seq(), expected_umi[i].as_slice());
            assert_eq!(rna_read.raw_umi_qual(), expected_umi_qual[i].as_slice());

            // Check that the R1 data is correct
            assert_eq!(rna_read.r1_range, RpRange::new(WhichRead::R2, 0, None));
            assert_eq!(rna_read.r1_seq(), expected_r1[i].as_slice());
            assert_eq!(rna_read.r1_qual(), expected_r1_qual[i].as_slice());

            // Check that the R2 data is correct
            assert_eq!(rna_read.r2_range, None);
            assert_eq!(rna_read.r2_seq(), None);
            assert_eq!(rna_read.r2_qual(), None);
        }
    }

    #[test]
    fn test_rna_chunk_subsample() {
        let chunk = RnaChunk {
            chemistry: get_vdj_chemistry(),
            gem_group: 1,
            read_chunks: ReadChunks{
                r1: "tests/rna_read/interleaved_2k.fastq".into(),
                r2: None,
                i1: None,
                i2: None,
            },
            read_group: "Blah".into(),
            reads_interleaved: true,
            subsample_rate: Some(0.2),
            read_lengths: FxHashMap::default(),
        };
        let processed_reads: usize = chunk.iter().unwrap().map(|_| 1).sum();
        println!("{}", processed_reads);
        // Expecting 400 reads since we start with 2000 reads
        assert!(processed_reads > 300); // < 1e-6 probability of this happening by chance
        assert!(processed_reads < 500); // < 1e-6 probability of this happening by chance
    }

    proptest! {
        #[test]
        fn prop_test_rna_chunk_trim(
            r1_length in 41usize..,
            trim_r1 in any::<bool>(),
            r2_length in any::<usize>(),
            trim_r2 in any::<bool>()
        ) { 
            // SC-VDJ chemistry
            {
                let expected_r1_length = if trim_r1 {
                    Some(min(150-41, r1_length-41))
                } else {
                    None
                };

                let expected_r2_length = if trim_r2 {
                    Some(min(150, r2_length))
                } else {
                    None
                };
                
                let mut chunk = RnaChunk {
                    chemistry: get_vdj_chemistry(),
                    gem_group: 1,
                    read_chunks: ReadChunks{
                        r1: "tests/rna_read/interleaved_2k.fastq.lz4".into(),
                        r2: None,
                        i1: None,
                        i2: None,
                    },
                    read_group: "Blah".into(),
                    reads_interleaved: true,
                    subsample_rate: None,
                    read_lengths: FxHashMap::default(),
                };
                if trim_r1 {
                    chunk.illumina_r1_trim_length(r1_length);
                }
                if trim_r2 {
                    chunk.illumina_r2_trim_length(r2_length);
                }

                for rna_read_result in chunk.iter().unwrap().take(100) {
                    let rna_read = rna_read_result.unwrap();
                    assert_eq!(rna_read.bc_range, RpRange::new(WhichRead::R1, 0, Some(16)));
                    assert_eq!(rna_read.umi_range, RpRange::new(WhichRead::R1, 16, Some(10)));
                    assert_eq!(rna_read.r1_range, RpRange::new(WhichRead::R1, 41, expected_r1_length));
                    assert_eq!(rna_read.r2_range, Some(RpRange::new(WhichRead::R2, 0, expected_r2_length)));
                }
            }
            // SC-VDJ R2 chemistry
            {
                let expected_r1_length = if trim_r2 {
                    Some(min(150, r2_length))
                } else {
                    None
                };

                let mut chunk = RnaChunk {
                    chemistry: get_vdj_r2_chemistry(),
                    gem_group: 1,
                    read_chunks: ReadChunks{
                        r1: "tests/rna_read/interleaved_2k.fastq.lz4".into(),
                        r2: None,
                        i1: None,
                        i2: None,
                    },
                    read_group: "Blah".into(),
                    reads_interleaved: true,
                    subsample_rate: None,
                    read_lengths: FxHashMap::default(),
                };
                if trim_r1 {
                    chunk.illumina_r1_trim_length(r1_length);
                }
                if trim_r2 {
                    chunk.illumina_r2_trim_length(r2_length);
                }

                for rna_read_result in chunk.iter().unwrap().take(100) {
                    let rna_read = rna_read_result.unwrap();
                    assert_eq!(rna_read.bc_range, RpRange::new(WhichRead::R1, 0, Some(16)));
                    assert_eq!(rna_read.umi_range, RpRange::new(WhichRead::R1, 16, Some(10)));
                    assert_eq!(rna_read.r1_range, RpRange::new(WhichRead::R2, 0, expected_r1_length));
                    assert_eq!(rna_read.r2_range, None);
                }
            }
        }
    }

    #[test]
    fn test_rna_read_adapter_trim() {
        // The VDJ adapters are listed in `tests/rna_read/vdj_adapters.json`. The reads in `tests/rna_read/interleaved_2k_insert.fastq` were
        // trimmed using cutadapt and saved in `tests/rna_read/interleaved_trimmed_2k.fastq` (See `tests/rna_read/run_cutadapt.sh`). This
        // test makes sure that the trimming that we produce exactly matches the cutadapt outputs. There are 2000 read pairs in total and 85
        // read pairs are trimmed in total.
        use adapter_trimmer::{Adapter, ReadAdapterCatalog};
        use read_pair_iter::ReadPairIter;

        let vdj_adapters: FxHashMap<WhichRead, Vec<Adapter>> = serde_json::from_reader(File::open("tests/rna_read/vdj_adapters.json").unwrap()).unwrap();
        let mut ad_catalog = ReadAdapterCatalog::new();
        for (read, adapters) in vdj_adapters.iter() {
            for adapter in adapters {
                ad_catalog.add_adapter(*read, adapter);
            }
        }

        let fastqs = InputFastqs {
            // This file was created by running the script `tests/rna_read/run_cutadapt.sh`
            // It is be the trimmed output fastq returned by cutadapt
            r1: "tests/rna_read/interleaved_trimmed_2k.fastq".into(),
            r2: None,
            i1: None,
            i2: None,
            r1_interleaved: true,
        };

        // SCVDJ Chemistry
        {
            let rp_iter = ReadPairIter::from_fastq_files(fastqs.clone()).unwrap();

            let chunk = RnaChunk {
                chemistry: get_vdj_chemistry(),
                gem_group: 1,
                read_chunks: ReadChunks{
                    r1: "tests/rna_read/interleaved_2k.fastq".into(),
                    r2: None,
                    i1: None,
                    i2: None,
                },
                read_group: "Blah".into(),
                reads_interleaved: true,
                subsample_rate: None,
                read_lengths: FxHashMap::default(),
            };

            let mut n_trimmed = 0;
            let mut adapter_counts: FxHashMap<String, usize> = FxHashMap::default();
            for (rna_read_result, rp_result) in chunk.iter().unwrap().zip(rp_iter) {
                let mut rna_read = rna_read_result.unwrap();
                let mut adapter_positions = rna_read.trim_adapters(&mut ad_catalog);
                for (k, _) in adapter_positions.drain() {
                    *adapter_counts.entry(k).or_insert(0) += 1;
                }
                let rp = rp_result.unwrap();
                assert_eq!(rna_read.r1_seq(), rp.get(WhichRead::R1, ReadPart::Seq).unwrap());
                assert_eq!(rna_read.r2_seq().unwrap(), rp.get(WhichRead::R2, ReadPart::Seq).unwrap());
                if rna_read.r1_seq().len() != 109 || rna_read.r2_seq().unwrap().len() != 150 {
                    n_trimmed += 1;
                }
            }
            println!("Trimmed {} sequences", n_trimmed);

            // The counts returned by cutadapt does not completely agree with this, because they count things
            // differently when multiple adapters are found for a read.
            assert_eq!(adapter_counts.get("R2_rc"), Some(&27));
            assert_eq!(adapter_counts.get("P7_rc"), Some(&8));
            assert_eq!(adapter_counts.get("polyA"), Some(&8));
            assert_eq!(adapter_counts.get("rt_primer_rc"), Some(&2));
            assert_eq!(adapter_counts.get("spacer"), None);
            assert_eq!(adapter_counts.get("spacer_rc"), Some(&66));
            assert_eq!(adapter_counts.get("R1_rc"), Some(&29));
            assert_eq!(adapter_counts.get("P5_rc"), Some(&14));
            assert_eq!(adapter_counts.get("polyT"), Some(&4));
            assert_eq!(adapter_counts.get("rt_primer"), Some(&1));
        }

        // SC-VDJ R2 chemistry
        {
            let rp_iter = ReadPairIter::from_fastq_files(fastqs.clone()).unwrap();

            let chunk = RnaChunk {
                chemistry: get_vdj_r2_chemistry(),
                gem_group: 1,
                read_chunks: ReadChunks{
                    r1: "tests/rna_read/interleaved_2k.fastq".into(),
                    r2: None,
                    i1: None,
                    i2: None,
                },
                read_group: "Blah".into(),
                reads_interleaved: true,
                subsample_rate: None,
                read_lengths: FxHashMap::default(),
            };

            for (rna_read_result, rp_result) in chunk.iter().unwrap().zip(rp_iter) {
                let mut rna_read = rna_read_result.unwrap();
                rna_read.trim_adapters(&mut ad_catalog);
                let rp = rp_result.unwrap();
                assert_eq!(rna_read.r1_seq(), rp.get(WhichRead::R2, ReadPart::Seq).unwrap());
                assert_eq!(rna_read.r2_seq(), None);
            }
        }
    }
}
