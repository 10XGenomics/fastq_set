

use Barcode;
use HasBarcode;
use shardio::*;
use failure::Error;
use metric::{Metric, SimpleHistogram, SerdeFormat};
use metric_utils::{ConfiguredReadMetric, ReadMetric, DefaultConfig};
use std::path::Path;
use serde::Serialize;
use barcode::BarcodeChecker;
use FastqProcessor;

// TODO: Set the parameters optimally?
pub const SEND_BUFFER_SZ: usize = 256;
pub const DISK_CHUNK_SZ: usize = 2048; // 2MB if ~1kB per record
pub const ITEM_BUFFER_SZ: usize = 2097152; // 2GB if above holds

pub struct BarcodeSortWorkflow<Processor, MetricType> 
where
    Processor: FastqProcessor,
    <Processor as FastqProcessor>::ReadType: 'static + HasBarcode + Send + Serialize,
    MetricType: ConfiguredReadMetric<ReadType=<Processor as FastqProcessor>::ReadType>,
{
    processor: Processor,
    sorter: BarcodeSorter<<Processor as FastqProcessor>::ReadType>,
    checker: BarcodeChecker,
    metrics: MetricType,
    config: <MetricType as ConfiguredReadMetric>::ConfigType,
}

impl<Processor, MetricType> BarcodeSortWorkflow<Processor, MetricType>
where
    Processor: FastqProcessor,
    <Processor as FastqProcessor>::ReadType: 'static + HasBarcode + Send + Serialize,
    MetricType: ReadMetric<ReadType=<Processor as FastqProcessor>::ReadType>,
{
    pub fn new<P: AsRef<Path>>(processor: Processor,
        valid_path: P, 
        invalid_path: P,
        whitelist_path: P) -> Result<Self, Error> {
        
        let sorter = BarcodeSorter::new(valid_path, invalid_path)?;
        let checker = BarcodeChecker::new(whitelist_path)?;
        Ok(BarcodeSortWorkflow {
            processor,
            sorter,
            checker,
            metrics: Metric::new(),
            config: DefaultConfig,
        })
    }
}

impl<Processor, MetricType> BarcodeSortWorkflow<Processor, MetricType>
where
    Processor: FastqProcessor,
    <Processor as FastqProcessor>::ReadType: 'static + HasBarcode + Send + Serialize,
    MetricType: ConfiguredReadMetric<ReadType=<Processor as FastqProcessor>::ReadType>,
{

    pub fn new_with_config<P: AsRef<Path>>(processor: Processor,
        valid_path: P,
        invalid_path: P,
        whitelist_path: P,
        config: <MetricType as ConfiguredReadMetric>::ConfigType) -> Result<Self, Error> {
        
        let sorter = BarcodeSorter::new(valid_path, invalid_path)?;
        let checker = BarcodeChecker::new(whitelist_path)?;
        Ok(BarcodeSortWorkflow {
            processor,
            sorter,
            checker,
            metrics: Metric::new(),
            config,
        })
    }

    pub fn execute_workflow(&mut self) -> Result<(), Error> {
        let mut nreads = 0;
        for read_result in self.processor.iter()? {
            nreads += 1;
            let mut read = read_result?;
            let mut bc = read.barcode().clone();
            self.checker.check(&mut bc);
            read.set_barcode(bc);
            self.metrics.update_with_config(&self.config, &read);
            self.sorter.process(read);
        }
        println!("Processed {} reads", nreads);
        self.sorter.finish()
    }

    pub fn write_results<P: AsRef<Path>>(&self,
        barcode_counts_path: P,
        metrics_path: P) {
        self.sorter.write_barcode_counts(barcode_counts_path);
        self.metrics.to_file(metrics_path, SerdeFormat::Binary);
    }
}

pub struct BarcodeOrder;
impl<T> SortKey<T, Barcode> for BarcodeOrder
where
    T: HasBarcode,
{
    fn sort_key(v: &T) -> &Barcode {
        v.barcode()
    }
}

struct BarcodeSorter<T>
where
    T: 'static + HasBarcode + Send + Serialize
{
    valid_writer: ShardWriter<T, Barcode, BarcodeOrder>,
    valid_sender: ShardSender<T>,
    invalid_writer: ShardWriter<T, Barcode, BarcodeOrder>,
    invalid_sender: ShardSender<T>,
    valid_bc_distribution: SimpleHistogram<Barcode>,
}

impl<T> BarcodeSorter<T>
where
    T: 'static + HasBarcode + Send + Serialize
{
    fn new<P: AsRef<Path>>(valid_path: P, invalid_path: P) -> Result<Self, Error> {
        let valid_writer: ShardWriter<T, Barcode, BarcodeOrder> = 
            ShardWriter::new(valid_path, SEND_BUFFER_SZ, DISK_CHUNK_SZ, ITEM_BUFFER_SZ)?;
        let valid_sender = valid_writer.get_sender();

        let invalid_writer: ShardWriter<T, Barcode, BarcodeOrder> = 
            ShardWriter::new(invalid_path, SEND_BUFFER_SZ, DISK_CHUNK_SZ, ITEM_BUFFER_SZ)?;
        let invalid_sender = invalid_writer.get_sender();

        Ok(BarcodeSorter {
            valid_writer,
            valid_sender,
            invalid_writer,
            invalid_sender,
            valid_bc_distribution: SimpleHistogram::new(),
        })
    }

    fn process(&mut self, read: T) {
        if read.barcode().is_valid() {
            self.valid_bc_distribution.observe(read.barcode().clone());
            self.valid_sender.send(read);
        } else {
            self.invalid_sender.send(read);
        }
    }

    fn finish(&mut self) -> Result<(), Error> {
        self.valid_sender.finished();
        self.invalid_sender.finished();
        self.valid_writer.finish()?;
        self.invalid_writer.finish()?;
        Ok(())
    }

    fn write_barcode_counts<P: AsRef<Path>>(&self, path: P) {
        self.valid_bc_distribution.to_file(path, SerdeFormat::Binary);
    }
}

pub fn write_merged_barcode_counts<P: AsRef<Path>, Q: AsRef<Path>>(shard_counts: &[P], merged_file: Q) {
    let bc_counts: SimpleHistogram<Barcode> = SimpleHistogram::from_files(shard_counts, SerdeFormat::Binary);
    bc_counts.to_file(merged_file, SerdeFormat::Binary);
}