#[macro_use]
extern crate criterion;
use criterion::Criterion;
use fastq_set::SSeq;

fn run_benchmark(c: &mut Criterion) {
    c.bench(
        "bench-byte-array",
        criterion::Benchmark::new("from-bytes", |b| {
            b.iter(|| SSeq::from_bytes(b"AGTCCTCTGCATTTTG"))
        })
        .with_function("from-bytes-unchecked", |b| {
            b.iter(|| SSeq::from_bytes_unchecked(b"AGTCCTCTGCATTTTG"))
        })
        .with_function("from-iter", |b| {
            b.iter(|| SSeq::from_iter(b"AGTCCTCTGCATTTTG"))
        }),
    );
}

criterion_group!(benches, run_benchmark);

criterion_main!(benches);
