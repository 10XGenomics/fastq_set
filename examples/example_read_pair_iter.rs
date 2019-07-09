extern crate fastq_10x;
use fastq_10x::read_pair_iter::ReadPairIter;

fn main() {
    let file = "tests/rna_read/interleaved_2k.fastq.lz4";
    let rp_iter = ReadPairIter::new(Some(file.to_string()), None, None, None, true).unwrap();
    let n_read_pairs = rp_iter.count();
    println!("{}", n_read_pairs);
}
