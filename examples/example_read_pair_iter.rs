extern crate fastq_10x;
use fastq_10x::read_pair_iter::ReadPairIter;

fn main() {
    let file = "/mnt/yard1/sreenath/tests/datasets/vdj_fastqs/micro_100k.fastq.lz4";
    let rp_iter = ReadPairIter::new(Some(file.to_string()), None, None, None, true).unwrap();
    let n_read_pairs = rp_iter.count();
    println!("{}", n_read_pairs);
}
