use serde_json;

use fastq_10x::read_pair_iter::InputFastqs;
use fastq_10x::rna_read::{ChemistryDef, RnaChunk};
use fastq_10x::FastqProcessor;
use std::fs::File;

fn main() {
    let chemistry: ChemistryDef = serde_json::from_reader(
        File::open("tests/rna_read/sc_vdj_chemistry.json").expect("Failed to open chemistry json"),
    )
    .unwrap();
    let file = "tests/rna_read/interleaved_2k.fastq";

    let chunk = RnaChunk::new(
        chemistry,
        1,
        InputFastqs {
            r1: file.to_string(),
            r2: None,
            i1: None,
            i2: None,
            r1_interleaved: true,
        },
        "my_group".into(),
    );
    let n_read_pairs = chunk.iter().unwrap().count();
    println!("{}", n_read_pairs);
}
