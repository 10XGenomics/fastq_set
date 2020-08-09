#![no_main]
#[macro_use] extern crate libfuzzer_sys;
extern crate fastq_set;

extern crate serde;
#[macro_use] extern crate serde_derive;
extern crate bincode;
use fastq_set::{Record, OwnedRecord};
use fastq_set::read_pair::{ReadPair, ReadPart, WhichRead};

#[derive(Debug, Serialize, Deserialize, Clone)]
struct Rec {
    hdr: Vec<u8>,
    seq: Vec<u8>,
}

fn map_rec_to_owned_record(rec: Rec) -> OwnedRecord {
    let alphabets = b"ACGT";
    let seq: Vec<u8> = rec.seq.into_iter().map(|x| {
        if x == b'N' {
            x 
        } else {
            alphabets[(x as usize) % alphabets.len()]
        }
    }).collect();

    let qual: Vec<u8> = seq.iter().map(|x| {
        if *x == b'N' {
            b'#' // Q score 2 
        } else {
            b'I' // Q40
        }
    }).collect();
    OwnedRecord {
        head: rec.hdr,
        seq,
        sep: None,
        qual,
    }
}

fn create_input(rr: [Option<Rec>; 4]) -> [Option<OwnedRecord>; 4] {
    let mut result: [Option<OwnedRecord>; 4] = [None, None, None, None];

    for (r_i, r_o) in rr.into_iter().zip(result.iter_mut()) {
        if r_i.is_some() {
            let r = r_i.clone().unwrap();
            *r_o = Some(map_rec_to_owned_record(r));
        }
    }
    result
}

use fastq::Parser;
fuzz_target!(|data: &[u8]| {
    if let Ok(rec) = bincode::deserialize::<[Option<Rec>; 4]>(data) {
        let input = create_input(rec);
        println!("Input : {:?}", input);

        // OwnedRecord does not implement Clone :/ So DIY
        let mut input_vec: Vec<Option<OwnedRecord>> = Vec::new();
        for r_opt in input.iter() {
            if let Some(r) = r_opt {
                input_vec.push(Some(OwnedRecord {
                    head: r.head().to_vec(),
                    seq: r.seq().to_vec(),
                    sep: None,
                    qual: r.qual.to_vec(),
                }));
            } else {
                input_vec.push(None);
            }
        }

        let read_pair = ReadPair::new(input);
        for (i, r) in input_vec.iter().enumerate() {
            let read = WhichRead::from(i);
            assert_eq!(read_pair.get(read, ReadPart::Header), r.as_ref().map(|x| x.head()));
            assert_eq!(read_pair.get(read, ReadPart::Seq), r.as_ref().map(|x| x.seq()));
            assert_eq!(read_pair.get(read, ReadPart::Qual), r.as_ref().map(|x| x.qual()));
        }
    }
});
