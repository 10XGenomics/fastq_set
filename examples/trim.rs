use fastq_10x::adapter_trimmer::{Adapter, AdapterTrimmer};

use failure::Error;

#[macro_use]
extern crate serde_derive;
use serde_json;

use std::fs::File;

use bio::io::fasta::{Reader, Writer};

use time::PreciseTime;

use std::env;

#[derive(Debug, Serialize, Deserialize)]
struct Input {
    adapter: Adapter,
    input_fasta: String,
    output_fasta: String,
}

fn main() -> Result<(), Error> {
    let start = PreciseTime::now();
    // Accepts a json which can populate the Input struct above
    let args: Vec<String> = env::args().collect();
    assert!(args.len() == 2);

    let input: Input = serde_json::from_reader(File::open(&args[1])?)?;

    let mut ntrimmed = 0;
    let mut total_reads = 0;

    let reader = Reader::from_file(&input.input_fasta)?;
    let mut writer = Writer::to_file(input.output_fasta)?;
    let mut trimmer = AdapterTrimmer::new(&input.adapter);

    for record in reader.records() {
        let record = record.unwrap();
        let seq = record.seq();

        let trimmed_seq;
        match trimmer.find(seq) {
            Some(result) => {
                ntrimmed += 1;
                trimmed_seq = &seq[result.retain_range];
            }
            None => {
                trimmed_seq = seq;
            }
        }
        writer.write(record.id(), record.desc(), trimmed_seq)?;
        total_reads += 1;
    }
    writer.flush()?;
    println!("Trimmed {}/{} reads", ntrimmed, total_reads);

    let end = PreciseTime::now();
    println!(
        "Elapsed time {:?} seconds for rust adapter trimmer.",
        start.to(end)
    );

    Ok(())
}
