use rust_htslib::bam::record::Record;
use {HasBarcode, Barcode, AlignableRead};
use bc_sort::BcSort;
use shardio::{ShardReaderSet, Range};
use serde::de::DeserializeOwned;
use itertools::Itertools;

// Assume PE for now.
pub trait DnaAligner {
    fn align(&self, sequence: &impl AlignableRead) -> Vec<(Record, Record)>;
}

struct PairAlign {
    prim_r1: Record,
    prim_r2: Record,

    supp_r1: Vec<Record>,
    supp_r2: Vec<Record>,
}

impl PairAlign {
    

}


pub trait ReadAlignMetrics<ReadType> {
    fn view(&self, reads: &[ReadType], alignments: &[Vec<(Record, Record)>]);
}

struct BcAlignProc<ReadType, A> {
    reader: ShardReaderSet<ReadType, Barcode, BcSort>,
    aligner: A,
}

impl<ReadType, A> BcAlignProc<ReadType, A> where
    ReadType: DeserializeOwned + HasBarcode + AlignableRead, A: DnaAligner {

    fn main(&self, chunk: Range<Barcode>) {
        let items = self.reader.iter_range(&chunk);

        for (bc, reads) in &items.group_by(|elt| elt.barcode()) {

            // Align the reads for this barcode
            let mut bc_alns = Vec::new();
            for r in reads {
                let alns = self.aligner.align(&r);
                bc_alns.push(alns);
            }

            // Mark duplicates (assumes duplicates have the same barcode)


            // Report read, alignment & barcode metrics.


            // Quantify 

            // Output

        }
    }
}