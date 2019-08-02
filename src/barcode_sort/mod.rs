use crate::barcode::BarcodeChecker;
use crate::Barcode;
use crate::FastqProcessor;
use crate::HasBarcode;
use failure::Error;
use metric::{Metric, SerdeFormat, SimpleHistogram};
use serde::Serialize;
use shardio::*;
use std;
use std::borrow::Cow;
use std::marker::PhantomData;
use std::path::Path;

// TODO: Set the parameters optimally?
pub const SEND_BUFFER_SZ: usize = 256;
pub const DISK_CHUNK_SZ: usize = 2048; // 2MB if ~1kB per record
pub const ITEM_BUFFER_SZ: usize = 2_097_152; // 2GB if above holds

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

pub struct BarcodeSortWorkflow<Processor, SortOrder>
where
    Processor: FastqProcessor,
    <Processor as FastqProcessor>::ReadType: 'static + HasBarcode + Send + Serialize,
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
        self.execute_workflow_with_visitor(max_iters, &mut v)
    }

    pub fn execute_workflow_with_visitor<
        Visitor: ReadVisitor<ReadType = <Processor as FastqProcessor>::ReadType>,
    >(
        &mut self,
        max_iters: Option<usize>,
        visitor: &mut Visitor,
    ) -> Result<(), Error> {
        let mut nreads = 0;
        for read_result in self
            .processor
            .iter()?
            .take(max_iters.unwrap_or(std::usize::MAX))
        {
            nreads += 1;
            let mut read = read_result?;
            let mut bc = *read.barcode();
            self.checker.check(&mut bc);
            read.set_barcode(bc);
            visitor.visit_read(&mut read)?;
            self.sorter.process(read);
        }
        println!("Processed {} reads", nreads);
        self.sorter.finish()
    }

    pub fn write_barcode_counts<P: AsRef<Path>>(&self, path: P) -> Result<(), Error> {
        self.sorter.write_barcode_counts(path)
    }

    pub fn num_valid_reads(&self) -> i64 {
        self.sorter.valid_reads
    }

    pub fn num_invalid_reads(&self) -> i64 {
        self.sorter.invalid_reads
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

struct BarcodeAwareSorter<T, Order>
where
    T: 'static + HasBarcode + Send + Serialize,
    Order: SortKey<T>,
{
    valid_writer: ShardWriter<T, Order>,
    valid_sender: ShardSender<T>,
    invalid_writer: ShardWriter<T, Order>,
    invalid_sender: ShardSender<T>,
    valid_bc_distribution: SimpleHistogram<Barcode>,
    valid_reads: i64,
    invalid_reads: i64,
}

impl<T, Order> BarcodeAwareSorter<T, Order>
where
    T: 'static + HasBarcode + Send + Serialize,
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
            valid_bc_distribution: SimpleHistogram::new(),
            valid_reads: 0,
            invalid_reads: 0,
        })
    }

    fn process(&mut self, read: T) {
        if read.barcode().is_valid() {
            self.valid_bc_distribution.observe(read.barcode());
            self.valid_sender.send(read);
            self.valid_reads += 1;
        } else {
            self.invalid_sender.send(read);
            self.invalid_reads += 1;
        }
    }

    fn finish(&mut self) -> Result<(), Error> {
        self.valid_sender.finished();
        self.invalid_sender.finished();
        self.valid_writer.finish()?;
        self.invalid_writer.finish()?;
        Ok(())
    }

    fn write_barcode_counts<P: AsRef<Path>>(&self, path: P) -> Result<(), Error> {
        self.valid_bc_distribution
            .to_file(path, SerdeFormat::Binary)
    }
}

pub fn write_merged_barcode_counts<P: AsRef<Path>, Q: AsRef<Path>>(
    shard_counts: &[P],
    merged_file: Q,
) -> Result<(), Error> {
    let bc_counts: SimpleHistogram<Barcode> =
        SimpleHistogram::from_files(shard_counts, SerdeFormat::Binary);
    bc_counts.to_file(merged_file, SerdeFormat::Binary)
}
