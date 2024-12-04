#[macro_use]
extern crate criterion;

use criterion::Criterion;
use fastq_set::read_pair_iter::ReadPairIter;
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

// fn rna_chunk_iter_count(chemistry: ChemistryDef, file: impl ToString) -> usize {
//     let chunk = RnaChunk::new(
//         chemistry,
//         1,
//         InputFastqs {
//             r1: file.to_string(),
//             r2: None,
//             i1: None,
//             i2: None,
//             r1_interleaved: true,
//         },
//         "my_group".into(),
//     );
//     chunk.iter().unwrap().count()
// }

// fn rna_read_trim(
//     chemistry: ChemistryDef,
//     fastq_file: impl ToString,
//     adapter_map: &HashMap<WhichRead, Vec<Adapter>>,
// ) -> usize {
//     let mut ad_catalog = ReadAdapterCatalog::from(adapter_map);
//     let chunk = RnaChunk::new(
//         chemistry,
//         1,
//         InputFastqs {
//             r1: fastq_file.to_string(),
//             r2: None,
//             i1: None,
//             i2: None,
//             r1_interleaved: true,
//         },
//         "my_group".into(),
//     );
//     let mut ntrimmed = 0;
//     for rna_read_result in chunk.iter().unwrap() {
//         let mut rna_read = rna_read_result.unwrap();
//         let adapter_positions = rna_read.trim_adapters(&mut ad_catalog);
//         if !adapter_positions.is_empty() {
//             ntrimmed += 1;
//         }
//     }
//     ntrimmed
// }

// fn cutadapt_trim(fastq_file: impl ToString) -> bool {
//     Command::new("sh")
//         .arg("time_cutadapt.sh")
//         .arg(fastq_file.to_string())
//         .status()
//         .is_ok()
// }

const INTERLEAVED_LZ4_FASTQ: &str = "tests/rna_read/interleaved_2k.fastq.lz4";
const INTERLEAVED_GZ_FASTQ: &str = "tests/rna_read/interleaved_2k.fastq.gz";
const INTERLEAVED_FASTQ: &str = "tests/rna_read/interleaved_2k.fastq";
// const INTERLEAVED_FASTQ_INSERT: &'static str = "tests/rna_read/interleaved_2k_insert.fastq";

fn run_fastq_lz4_benchmark(c: &mut Criterion) {
    c.bench_function("bench-lz4-wc", |b| {
        b.iter(
            || assert_eq!(lz4_count(INTERLEAVED_LZ4_FASTQ), 16000), // 16000 lines
        );
    });
    c.bench_function("bench-lz4-read-pair-iter-count", |b| {
        b.iter(
            || assert_eq!(read_pair_count(INTERLEAVED_LZ4_FASTQ), 2000), // 2000 read pairs
        );
    });
    // let chemistry: ChemistryDef =
    //     serde_json::from_reader(File::open("tests/rna_read/sc_vdj_chemistry.json").unwrap())
    //         .unwrap();
    // c.bench_function("bench-lz4-rna-chunk-iter-count", move |b| {
    //     b.iter(
    //         || {
    //             assert_eq!(
    //                 rna_chunk_iter_count(chemistry.clone(), INTERLEAVED_LZ4_FASTQ),
    //                 2000
    //             )
    //         }, // 2000 read pairs
    //     )
    // });
}

fn run_fastq_gz_benchmark(c: &mut Criterion) {
    c.bench_function("bench-gz-wc", |b| {
        b.iter(
            || assert_eq!(gz_count(INTERLEAVED_GZ_FASTQ), 16000), // 16000 lines
        );
    });
    c.bench_function("bench-gz-read-pair-iter-count", |b| {
        b.iter(
            || assert_eq!(read_pair_count(INTERLEAVED_GZ_FASTQ), 2000), // 2000 read pairs
        );
    });
    // let chemistry: ChemistryDef =
    //     serde_json::from_reader(File::open("tests/rna_read/sc_vdj_chemistry.json").unwrap())
    //         .unwrap();
    // c.bench_function("bench-gz-rna-chunk-iter-count", move |b| {
    //     b.iter(
    //         || {
    //             assert_eq!(
    //                 rna_chunk_iter_count(chemistry.clone(), INTERLEAVED_GZ_FASTQ),
    //                 2000
    //             )
    //         }, // 2000 read pairs
    //     )
    // });
}

fn run_fastq_benchmark(c: &mut Criterion) {
    c.bench_function("bench-fastq-wc", |b| {
        b.iter(
            || assert_eq!(simple_count(INTERLEAVED_FASTQ), 16000), // 16000 lines
        );
    });
    c.bench_function("bench-fastq-read-pair-iter-count", |b| {
        b.iter(
            || assert_eq!(read_pair_count(INTERLEAVED_FASTQ), 2000), // 2000 read pairs
        );
    });
    // let chemistry: ChemistryDef =
    //     serde_json::from_reader(File::open("tests/rna_read/sc_vdj_chemistry.json").unwrap())
    //         .unwrap();
    // c.bench_function("bench-fastq-rna-chunk-iter-count", move |b| {
    //     b.iter(
    //         || {
    //             assert_eq!(
    //                 rna_chunk_iter_count(chemistry.clone(), INTERLEAVED_FASTQ),
    //                 2000
    //             )
    //         }, // 2000 read pairs
    //     )
    // });
}

// fn run_trim_benchmark(c: &mut Criterion) {
//     let chemistry: ChemistryDef =
//         serde_json::from_reader(File::open("tests/rna_read/sc_vdj_chemistry.json").unwrap())
//             .unwrap();
//     let vdj_adapters: HashMap<WhichRead, Vec<Adapter>> =
//         serde_json::from_reader(File::open("tests/rna_read/vdj_adapters.json").unwrap()).unwrap();

//     c.bench(
//         "bench-read-trim",
//         Benchmark::new("rust-utils", move |b| {
//             b.iter(|| rna_read_trim(chemistry.clone(), INTERLEAVED_FASTQ, &vdj_adapters))
//         })
//         .with_function("cutadapt", |b| {
//             b.iter(|| assert!(cutadapt_trim(INTERLEAVED_FASTQ_INSERT)))
//         })
//         .sample_size(10)
//         .throughput(Throughput::Elements(2000)),
//     );
// }

fn sseq_serde_bincode(c: &mut Criterion) {
    c.bench_function("bench-sseq-serde", |b| {
        let seq = b"AGCTAGTCAGTCAGTA";
        let mut sseqs = Vec::new();
        for _i in 0..10_000 {
            let s = fastq_set::sseq::SSeq::from_bytes(seq);
            sseqs.push(s);
        }

        b.iter(|| {
            let mut b = Vec::new();
            bincode::serialize_into(&mut b, &sseqs).unwrap();
            assert!(!b.is_empty());
        });
    });
}

criterion_group!(
    benches,
    run_fastq_lz4_benchmark,
    run_fastq_benchmark,
    run_fastq_gz_benchmark,
    // run_trim_benchmark,
    sseq_serde_bincode,
);

criterion_main!(benches);
