# fastq_set

Tools for working with groups FASTQ of files. 
Major functionality includes:
 * Find groups FASTQs (R1/R2/I1/I2) following Illumina filename conventions
 * Parsing flowcell information from Illumina FASTQ headers
 * High-speed FASTQ I/O (via the `fastq` crate), with careful validation of FASTQ correctness and good error message.
 * Containers for FASTQ read-pairs (along with index reads), providing access to 'technical' read components like cell barcode and
 * UMI sequences.
 * Flexible read trimming inspired by `cutadapt`