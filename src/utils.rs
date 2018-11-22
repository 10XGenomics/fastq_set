// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

//! Utility methods.

use std::boxed::Box;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter};
use std::path::Path;

use std::fmt::Debug;

use bincode;
use bincode::{deserialize_from, serialize_into};
use serde::de::DeserializeOwned;
use serde::Serialize;

use failure::Error;
use flate2::read::MultiGzDecoder;
use lz4;

const GZ_BUF_SIZE: usize = 1 << 22;

/// Open a (possibly gzipped) file into a BufReader.
pub fn open_with_gz<P: AsRef<Path>>(p: P) -> Result<Box<BufRead>, Error> {
    let r = File::open(p.as_ref())?;

    let ext = p.as_ref().extension().unwrap();

    if ext == "gz" {
        let gz = MultiGzDecoder::new(r);
        let buf_reader = BufReader::with_capacity(GZ_BUF_SIZE, gz);
        Ok(Box::new(buf_reader))
    } else if ext == "lz4" {
        let lz = lz4::Decoder::new(r)?;
        let buf_reader = BufReader::with_capacity(GZ_BUF_SIZE, lz);
        Ok(Box::new(buf_reader))
    } else {
        let buf_reader = BufReader::with_capacity(32 * 1024, r);
        Ok(Box::new(buf_reader))
    }
}

/// Serialize object `obj` of type `T` to the file `filename`
pub fn write_obj<T: Serialize, P: AsRef<Path> + Debug>(
    obj: &T,
    filename: P,
) -> Result<(), bincode::Error> {
    let f = match File::create(&filename) {
        Err(err) => panic!("couldn't create file {:?}: {}", filename, err),
        Ok(f) => f,
    };
    let mut writer = BufWriter::new(f);
    serialize_into(&mut writer, &obj)
}

/// Serialize an object  of type `T` from the file `filename`
pub fn read_obj<T: DeserializeOwned, P: AsRef<Path> + Debug>(
    filename: P,
) -> Result<T, bincode::Error> {
    let f = match File::open(&filename) {
        Err(err) => panic!("couldn't open file {:?}: {}", filename, err),
        Ok(f) => f,
    };
    let mut reader = BufReader::new(f);
    deserialize_from(&mut reader)
}
