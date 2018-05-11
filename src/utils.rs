//
// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
//

//! Various useful methods (could be cleaned up)

use std::fs::File;
use std::path::{Path};
use std::io::{BufReader, BufWriter};
use std::io::BufRead;

use std::fmt::Debug;

use bincode;
use serde::{Serialize};
use serde::de::DeserializeOwned;
use bincode::{serialize_into, deserialize_from};



pub fn write_obj<T: Serialize, P: AsRef<Path> + Debug>(g: &T, filename: P) -> Result<(), bincode::Error> {
    let f = match File::create(&filename) {
        Err(err) => panic!("couldn't create file {:?}: {}", filename, err),
        Ok(f) => f,
    };
    let mut writer = BufWriter::new(f);
    serialize_into(&mut writer, &g)
}

pub fn read_obj<T: DeserializeOwned, P: AsRef<Path> + Debug>(filename: P) -> Result<T, bincode::Error> {
    let f = match File::open(&filename) {
        Err(err) => panic!("couldn't open file {:?}: {}", filename, err),
        Ok(f) => f,
    };
    let mut reader = BufReader::new(f);
    deserialize_from(&mut reader)
}

