use csv;
use csv::WriterBuilder;
use failure::Error;
use rayon;
use rayon::prelude::*;
use serde::Serialize;
use std::collections::HashMap;
use std::path::Path;

use fastq_10x::filenames::bcl_processor::{find_flowcell_fastqs, BclProcessorFile};
use fastq_10x::read_pair::{ReadPart::Seq, WhichRead::I2};
use fastq_10x::read_pair_iter::InputFastqs;
use fastq_10x::sseq::SSeq;

pub fn count_dual_indexes(fastqs: InputFastqs) -> Result<HashMap<SSeq, u32>, Error> {
    let read_pair_iter = fastq_10x::read_pair_iter::ReadPairIter::from_fastq_files(&fastqs)?;
    let mut counts = HashMap::new();
    let mut n = 0;

    for r in read_pair_iter {
        let r = r?;

        let i2 = SSeq::new(r.get(I2, Seq).unwrap());
        let c = counts.entry(i2).or_insert(0);
        *c += 1;
        n += 1;
    }

    println!("read {} records from {:?}", n, fastqs.i1);
    Ok(counts)
}

pub fn count_dual_indexes_flowcell(
    p: impl AsRef<Path>,
) -> Result<Vec<(BclProcessorFile, HashMap<SSeq, u32>)>, Error> {
    let fc_data = find_flowcell_fastqs(p)?;
    let mut results = Vec::new();

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(8)
        .build()
        .unwrap();

    pool.install(|| {
        fc_data
            .into_par_iter()
            .map(|(info, fastqs)| count_dual_indexes(fastqs).map(|f| (info, f)))
            .collect_into_vec(&mut results);
    });

    results.into_iter().collect()
}

#[derive(Serialize)]
struct CountStruct<'a, 'b> {
    lane: usize,
    chunk: usize,
    i1: &'a str,
    i2: &'b str,
    n: u32,
}

pub fn write_counts(
    file: impl AsRef<Path>,
    data: Vec<(BclProcessorFile, HashMap<SSeq, u32>)>,
) -> Result<(), Error> {
    let writer = std::io::BufWriter::new(std::fs::File::create(file)?);
    let mut wtr = WriterBuilder::new().from_writer(writer);

    for (bcl, counts) in data {
        for (i2, count) in counts {
            let cs = CountStruct {
                lane: bcl.lane,
                chunk: bcl.chunk,
                i1: &bcl.si,
                i2: std::str::from_utf8(i2.seq()).unwrap(),
                n: count,
            };

            wtr.serialize(&cs)?;
        }
    }

    Ok(())
}

fn main() -> Result<(), Error> {
    let mut args = std::env::args();

    args.next();

    let p = args.next().unwrap();

    println!("Processing: {}", p);

    let counts = count_dual_indexes_flowcell(p)?;
    write_counts("counts.csv", counts)?;

    Ok(())
}
