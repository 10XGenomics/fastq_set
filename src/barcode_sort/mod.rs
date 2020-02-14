use crate::barcode::BarcodeChecker;
use crate::Barcode;
use crate::FastqProcessor;
use crate::HasBarcode;
use failure::Error;
use metric::{Metric, SerdeFormat, SimpleHistogram, TxHashMap};
use serde::{Deserialize, Serialize};
use shardio::*;
use std;
use std::borrow::Cow;
use std::hash::Hash;
use std::marker::PhantomData;
use std::path::Path;

// TODO: Set the parameters optimally?
pub const SEND_BUFFER_SZ: usize = 256;
pub const DISK_CHUNK_SZ: usize = 8_192; // 8MiB if ~1kB per record
pub const ITEM_BUFFER_SZ: usize = 1_048_576; // 1GiB if above holds

pub trait ReadVisitor {
    type ReadType;
    fn visit_read(&mut self, read: &mut Self::ReadType) -> Result<(), Error>;
}

struct DefaultVisitor<T> {
    phantom: PhantomData<T>,
}
impl<T> ReadVisitor for DefaultVisitor<T> {
    type ReadType = T;
    fn visit_read(&mut self, _: &mut Self::ReadType) -> Result<(), Error> {
        Ok(())
    }
}

pub struct BarcodeSortWorkflow<Processor, SortOrder = BarcodeOrder>
where
    Processor: FastqProcessor,
    <Processor as FastqProcessor>::ReadType: 'static + HasBarcode + Send + Serialize,
    <<Processor as FastqProcessor>::ReadType as HasBarcode>::LibraryType:
        Eq + Hash + Clone + Serialize + for<'de> Deserialize<'de>,
    SortOrder: SortKey<<Processor as FastqProcessor>::ReadType>,
    <SortOrder as SortKey<<Processor as FastqProcessor>::ReadType>>::Key:
        'static + Send + Ord + Serialize + Clone,
{
    processor: Processor,
    sorter: BarcodeAwareSorter<<Processor as FastqProcessor>::ReadType, SortOrder>,
    checker: BarcodeChecker,
}

impl<Processor, SortOrder> BarcodeSortWorkflow<Processor, SortOrder>
where
    Processor: FastqProcessor,
    <Processor as FastqProcessor>::ReadType: 'static + HasBarcode + Send + Serialize,
    <<Processor as FastqProcessor>::ReadType as HasBarcode>::LibraryType:
        Eq + Hash + Clone + Serialize + for<'de> Deserialize<'de>,
    SortOrder: SortKey<<Processor as FastqProcessor>::ReadType>,
    <SortOrder as SortKey<<Processor as FastqProcessor>::ReadType>>::Key:
        'static + Send + Ord + Serialize + Clone,
{
    pub fn new(
        processor: Processor,
        valid_path: impl AsRef<Path>,
        invalid_path: impl AsRef<Path>,
        whitelist_path: impl AsRef<Path>,
    ) -> Result<Self, Error> {
        let sorter = BarcodeAwareSorter::new(valid_path, invalid_path)?;
        let checker = BarcodeChecker::new(whitelist_path)?;
        Ok(BarcodeSortWorkflow {
            processor,
            sorter,
            checker,
        })
    }

    pub fn execute_workflow(&mut self, max_iters: Option<usize>) -> Result<(), Error> {
        let mut v = DefaultVisitor {
            phantom: PhantomData,
        };
        // FIXME: hoist to caller?
        let f = |_: &<Processor as FastqProcessor>::ReadType| false;
        self.execute_workflow_with_visitor(max_iters, &mut v, f)
    }

    pub fn execute_workflow_with_visitor<V, F>(
        &mut self,
        max_iters: Option<usize>,
        visitor: &mut V,
        translate_barcode: F,
    ) -> Result<(), Error>
    where
        V: ReadVisitor<ReadType = <Processor as FastqProcessor>::ReadType>,
        F: Fn(&<Processor as FastqProcessor>::ReadType) -> bool,
    {
        let mut nreads = 0;
        for read_result in self
            .processor
            .iter()?
            .take(max_iters.unwrap_or(std::usize::MAX))
        {
            nreads += 1;
            let mut read = read_result?;
            let mut bc = *read.barcode();
            self.checker.check(&mut bc, translate_barcode(&read));
            read.set_barcode(bc);
            visitor.visit_read(&mut read)?;
            self.sorter.process(read)?;
        }
        println!("Processed {} reads", nreads);
        self.sorter.finish()
    }

    pub fn write_barcode_counts<P: AsRef<Path>>(&self, path: P) -> Result<(), Error> {
        self.sorter.write_barcode_counts(path)
    }

    pub fn num_valid_items(
        &self,
    ) -> &TxHashMap<<<Processor as FastqProcessor>::ReadType as HasBarcode>::LibraryType, i64> {
        &self.sorter.valid_items
    }

    pub fn num_invalid_items(
        &self,
    ) -> &TxHashMap<<<Processor as FastqProcessor>::ReadType as HasBarcode>::LibraryType, i64> {
        &self.sorter.invalid_items
    }
}

pub struct BarcodeOrder;
impl<T> SortKey<T> for BarcodeOrder
where
    T: HasBarcode,
{
    type Key = Barcode;
    fn sort_key(v: &T) -> Cow<'_, Barcode> {
        Cow::Borrowed(v.barcode())
    }
}

struct BarcodeAwareSorter<T, Order = BarcodeOrder>
where
    T: 'static + HasBarcode + Send + Serialize,
    <T as HasBarcode>::LibraryType: Eq + Hash + Clone + Serialize + for<'de> Deserialize<'de>,
    Order: SortKey<T>,
    Order::Key: 'static + Send + Serialize,
{
    // note: sender must appear before writers, because senders
    // must be dropped before writers.
    valid_sender: ShardSender<T, Order>,
    valid_writer: ShardWriter<T, Order>,
    invalid_sender: ShardSender<T, Order>,
    invalid_writer: ShardWriter<T, Order>,
    valid_bc_distribution: TxHashMap<<T as HasBarcode>::LibraryType, SimpleHistogram<Barcode>>,
    valid_items: TxHashMap<<T as HasBarcode>::LibraryType, i64>,
    invalid_items: TxHashMap<<T as HasBarcode>::LibraryType, i64>,
}

impl<T, Order> BarcodeAwareSorter<T, Order>
where
    T: 'static + HasBarcode + Send + Serialize,
    <T as HasBarcode>::LibraryType: Eq + Hash + Clone + Serialize + for<'de> Deserialize<'de>,
    Order: SortKey<T>,
    <Order as SortKey<T>>::Key: 'static + Send + Ord + Serialize + Clone,
{
    fn new(valid_path: impl AsRef<Path>, invalid_path: impl AsRef<Path>) -> Result<Self, Error> {
        let valid_writer: ShardWriter<T, Order> =
            ShardWriter::new(valid_path, SEND_BUFFER_SZ, DISK_CHUNK_SZ, ITEM_BUFFER_SZ)?;
        let valid_sender = valid_writer.get_sender();

        let invalid_writer: ShardWriter<T, Order> =
            ShardWriter::new(invalid_path, SEND_BUFFER_SZ, DISK_CHUNK_SZ, ITEM_BUFFER_SZ)?;
        let invalid_sender = invalid_writer.get_sender();

        Ok(BarcodeAwareSorter {
            valid_writer,
            valid_sender,
            invalid_writer,
            invalid_sender,
            valid_bc_distribution: TxHashMap::default(),
            valid_items: TxHashMap::default(),
            invalid_items: TxHashMap::default(),
        })
    }

    fn process(&mut self, read: T) -> Result<(), Error> {
        if read.barcode().is_valid() {
            self.valid_bc_distribution
                .entry(read.library_type())
                .or_insert(SimpleHistogram::new())
                .observe(read.barcode());
            *self.valid_items.entry(read.library_type()).or_insert(0) += 1;
            self.valid_sender.send(read)?;
        } else {
            *self.invalid_items.entry(read.library_type()).or_insert(0) += 1;
            self.invalid_sender.send(read)?;
        }
        Ok(())
    }

    fn finish(&mut self) -> Result<(), Error> {
        self.valid_sender.finished()?;
        self.invalid_sender.finished()?;
        self.valid_writer.finish()?;
        self.invalid_writer.finish()?;
        Ok(())
    }

    fn write_barcode_counts<P: AsRef<Path>>(&self, path: P) -> Result<(), Error> {
        self.valid_bc_distribution
            .to_file(path, SerdeFormat::Binary)
    }
}

pub fn write_merged_barcode_counts<P, Q, R>(shard_counts: &[P], merged_file: Q) -> Result<(), Error>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
    R: Eq + Hash + Clone + Serialize + for<'de> Deserialize<'de>,
{
    let bc_counts: TxHashMap<R, SimpleHistogram<Barcode>> = {
        let mut res = TxHashMap::default();
        for x in shard_counts {
            let y = std::io::BufReader::new(std::fs::File::open(x)?);
            let mut z: TxHashMap<R, SimpleHistogram<Barcode>> = bincode::deserialize_from(y)?;
            for (typ, hist) in z.drain() {
                res.entry(typ).or_insert(SimpleHistogram::new()).merge(hist)
            }
        }
        res
    };
    let handle = std::io::BufWriter::new(std::fs::File::create(merged_file)?);
    bincode::serialize_into(handle, &bc_counts)?;
    Ok(())
}
