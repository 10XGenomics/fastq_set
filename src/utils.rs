// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

//! Utility methods.

use std::any::{Any, TypeId};
use std::boxed::Box;
use std::collections::hash_map::DefaultHasher;
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::{BufRead, BufReader};
use std::path::Path;

use std::fmt::Debug;

use bincode::{deserialize_from, serialize_into};
use serde::de::DeserializeOwned;
use serde::Serialize;

use failure::{format_err, Error};
use flate2::read::MultiGzDecoder;
use lz4;

const GZ_BUF_SIZE: usize = 1 << 22;

/// Open a (possibly gzipped) file into a BufReader.
pub fn open_with_gz<P: AsRef<Path>>(p: P) -> Result<Box<dyn BufRead>, Error> {
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
pub fn write_obj<T: Any + Serialize, P: AsRef<Path> + Debug>(
    obj: &T,
    filename: P,
) -> Result<(), Error> {
    let f = File::create(&filename)?;
    let mut writer = lz4::EncoderBuilder::new().build(f)?;

    let typeid = TypeId::of::<T>();
    let mut hasher = DefaultHasher::new();
    typeid.hash(&mut hasher);
    let type_hash = hasher.finish();

    serialize_into(&mut writer, &type_hash)?;
    serialize_into(&mut writer, &obj)?;

    let (_, result) = writer.finish();
    result?;
    Ok(())
}

/// Serialize an object  of type `T` from the file `filename`
pub fn read_obj<T: Any + DeserializeOwned, P: AsRef<Path> + Debug>(
    filename: P,
) -> Result<T, Error> {
    let f = File::open(&filename)?;
    let lz4_reader = lz4::Decoder::new(f)?;
    let mut reader = BufReader::new(lz4_reader);

    let typeid = TypeId::of::<T>();
    let mut hasher = DefaultHasher::new();
    typeid.hash(&mut hasher);
    let type_hash = hasher.finish();

    let file_type_hash: u64 = deserialize_from(&mut reader)?;

    if type_hash != file_type_hash {
        let err = format_err!(
            "data type in file '{:?}' does not match expected type",
            filename.as_ref()
        );
        return Err(err);
    }

    let res = deserialize_from(&mut reader)?;
    Ok(res)
}

#[cfg(test)]
mod test {
    use crate::utils::*;
    use serde::{Deserialize, Serialize};

    #[test]
    fn test_write_obj_read_obj() -> Result<(), Error> {
        let fn1 = "test1.bin";
        let fn2 = "test2.bin";

        write_obj(&vec![1u8, 2u8, 3u8], fn1).unwrap();
        write_obj(&vec![1u64, 2u64, 3u64], fn2).unwrap();

        let _read1: Vec<u8> = read_obj(fn1)?;
        let _read2: Vec<u64> = read_obj(fn2)?;

        let read1_wrong: Result<Vec<u64>, _> = read_obj(fn1);
        assert!(read1_wrong.is_err());

        let read2_wrong: Result<Vec<u8>, _> = read_obj(fn2);
        assert!(read2_wrong.is_err());

        Ok(())
    }

    #[derive(Debug, Deserialize, Serialize, PartialEq, Eq)]
    struct T1 {
        arr1: Vec<(u16, u64, String)>,
        arr2: Vec<(String, u32, u8)>,
    }

    #[test]
    fn test_big_roundtrip() -> Result<(), Error> {
        let fn1 = "test1.bin";

        let mut arr1 = Vec::new();
        let mut arr2 = Vec::new();

        for i in 0..30000 {
            let v1 = (i % 5000) as u16;
            let v2 = i * 100;
            let v3 = "werqwer".to_string();
            arr1.push((v1, v2, v3));
        }

        for i in 0..50000 {
            let v1 = "eqwrdv".to_string();
            let v2 = (i * 100) as u32;
            let v3 = (i % 100) as u8;
            arr2.push((v1, v2, v3));
        }

        let obj = T1 { arr1, arr2 };

        write_obj(&obj, fn1)?;
        let read_obj: T1 = read_obj(fn1)?;

        assert_eq!(obj, read_obj);
        Ok(())
    }

}
