#[macro_use]
extern crate criterion;
extern crate fastq_10x;

use criterion::Criterion;
use fastq_10x::read_pair_iter::{InputFastqs, ReadPairIter};
use fastq_10x::rna_read::{ChemistryDef, RnaChunk};
use fastq_10x::FastqProcessor;
use std::fs::File;
use std::process::Command;

fn simple_count(file: impl ToString) -> usize {
    let stdout = String::from_utf8(
        Command::new("wc")
            .arg("-l")
            .arg(file.to_string())
            .output()
            .unwrap()
            .stdout
            .split(|x| *x == b' ')
            .next()
            .unwrap()
            .to_vec(),
    )
    .unwrap();
    stdout.trim().parse::<usize>().unwrap()
}

fn lz4_count(file: impl ToString) -> usize {
    let stdout = String::from_utf8(
        Command::new("sh")
            .arg("lz4_wc.sh")
            .arg(file.to_string())
            .output()
            .unwrap()
            .stdout,
    )
    .unwrap();
    stdout.trim().parse::<usize>().unwrap()
}

fn gz_count(file: impl ToString) -> usize {
    let stdout = String::from_utf8(
        Command::new("sh")
            .arg("gz_wc.sh")
            .arg(file.to_string())
            .output()
            .unwrap()
            .stdout,
    )
    .unwrap();
    stdout.trim().parse::<usize>().unwrap()
}

fn read_pair_count(file: impl ToString) -> usize {
    let rp_iter = ReadPairIter::new(Some(file.to_string()), None, None, None, true).unwrap();
    rp_iter.count()
}

fn rna_chunk_iter_count(chemistry: ChemistryDef, file: impl ToString) -> usize {
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
    chunk.iter().unwrap().count()
}

const INTERLEAVED_LZ4_FASTQ: &'static str = "tests/rna_read/interleaved_2k.fastq.lz4";
const INTERLEAVED_GZ_FASTQ: &'static str = "tests/rna_read/interleaved_2k.fastq.gz";
const INTERLEAVED_FASTQ: &'static str = "tests/rna_read/interleaved_2k.fastq";

fn run_fastq_lz4_benchmark(c: &mut Criterion) {
    c.bench_function("bench-lz4-wc", |b| {
        b.iter(
            || assert_eq!(lz4_count(INTERLEAVED_LZ4_FASTQ), 16000), // 16000 lines
        )
    });
    c.bench_function("bench-lz4-read-pair-iter-count", |b| {
        b.iter(
            || assert_eq!(read_pair_count(INTERLEAVED_LZ4_FASTQ), 2000), // 2000 read pairs
        )
    });
    let chemistry: ChemistryDef =
        serde_json::from_reader(File::open("tests/rna_read/sc_vdj_chemistry.json").unwrap())
            .unwrap();
    c.bench_function("bench-lz4-rna-chunk-iter-count", move |b| {
        b.iter(
            || {
                assert_eq!(
                    rna_chunk_iter_count(chemistry.clone(), INTERLEAVED_LZ4_FASTQ),
                    2000
                )
            }, // 2000 read pairs
        )
    });
}

fn run_fastq_gz_benchmark(c: &mut Criterion) {
    c.bench_function("bench-gz-wc", |b| {
        b.iter(
            || assert_eq!(gz_count(INTERLEAVED_GZ_FASTQ), 16000), // 16000 lines
        )
    });
    c.bench_function("bench-gz-read-pair-iter-count", |b| {
        b.iter(
            || assert_eq!(read_pair_count(INTERLEAVED_GZ_FASTQ), 2000), // 2000 read pairs
        )
    });
    let chemistry: ChemistryDef =
        serde_json::from_reader(File::open("tests/rna_read/sc_vdj_chemistry.json").unwrap())
            .unwrap();
    c.bench_function("bench-gz-rna-chunk-iter-count", move |b| {
        b.iter(
            || {
                assert_eq!(
                    rna_chunk_iter_count(chemistry.clone(), INTERLEAVED_GZ_FASTQ),
                    2000
                )
            }, // 2000 read pairs
        )
    });
}

fn run_fastq_benchmark(c: &mut Criterion) {
    c.bench_function("bench-fastq-wc", |b| {
        b.iter(
            || assert_eq!(simple_count(INTERLEAVED_FASTQ), 16000), // 16000 lines
        )
    });
    c.bench_function("bench-fastq-read-pair-iter-count", |b| {
        b.iter(
            || assert_eq!(read_pair_count(INTERLEAVED_FASTQ), 2000), // 2000 read pairs
        )
    });
    let chemistry: ChemistryDef =
        serde_json::from_reader(File::open("tests/rna_read/sc_vdj_chemistry.json").unwrap())
            .unwrap();
    c.bench_function("bench-fastq-rna-chunk-iter-count", move |b| {
        b.iter(
            || {
                assert_eq!(
                    rna_chunk_iter_count(chemistry.clone(), INTERLEAVED_FASTQ),
                    2000
                )
            }, // 2000 read pairs
        )
    });
}

criterion_group!(
    benches,
    run_fastq_lz4_benchmark,
    run_fastq_benchmark,
    run_fastq_gz_benchmark
);

criterion_main!(benches);
