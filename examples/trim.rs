use fastq_set::adapter_trimmer::{Adapter, AdapterTrimmer};

use anyhow::Error;

#[macro_use]
extern crate serde_derive;

use std::fs::File;

use bio::io::fasta::{Reader, Writer};

use std::time::Instant;

use std::env;

#[derive(Debug, Serialize, Deserialize)]
struct Input {
    adapter: Adapter,
    input_fasta: String,
    output_fasta: String,
}

fn main() -> Result<(), Error> {
    let start = Instant::now();
    // Accepts a json which can populate the Input struct above
    let args: Vec<String> = env::args().collect();
    assert!(args.len() == 2);

    let input: Input = serde_json::from_reader(File::open(&args[1])?)?;

    let mut ntrimmed = 0;
    let mut total_reads = 0;

    let reader = Reader::from_file(&input.input_fasta).unwrap();
    let mut writer = Writer::to_file(input.output_fasta)?;
    let mut trimmer = AdapterTrimmer::new(&input.adapter);

    for record in reader.records() {
        let record = record.unwrap();
        let seq = record.seq();

        let trimmed_seq = match trimmer.find(seq) {
            Some(result) => {
                ntrimmed += 1;
                &seq[result.retain_range]
            }
            None => seq,
        };
        writer.write(record.id(), record.desc(), trimmed_seq)?;
        total_reads += 1;
    }
    writer.flush()?;
    println!("Trimmed {ntrimmed}/{total_reads} reads");

    let end = Instant::now();
    println!(
        "Elapsed time {:?} seconds for rust adapter trimmer.",
        (end - start).as_secs_f32()
    );

    Ok(())
}
