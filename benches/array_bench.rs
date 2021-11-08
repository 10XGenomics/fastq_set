#[macro_use]
extern crate criterion;
use criterion::Criterion;
use fastq_set::SSeq;

fn run_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("bench-byte-array");
    group.bench_function("from-bytes", |b| {
        b.iter(|| SSeq::from_bytes(b"AGTCCTCTGCATTTTG"))
    });
    group.bench_function("from-bytes-unchecked", |b| {
        b.iter(|| SSeq::from_bytes_unchecked(b"AGTCCTCTGCATTTTG"))
    });
    group.bench_function("from-iter", |b| {
        b.iter(|| SSeq::from_iter(b"AGTCCTCTGCATTTTG"))
    });
    group.finish();
}

criterion_group!(benches, run_benchmark);

criterion_main!(benches);
