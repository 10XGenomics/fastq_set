# Streamlined Front-End Project

## Why do this?
- Reduce pipeline resource consumption: improve customer experience (reliability & cost), improve developer productivity
- Engineering load: better code reuse across pipelines, avoid continued copy/paste of front-end code
- Platform for future experimentation: new aligners, etc.
- Others??


## High-level Goals
- Reusable toolkit applicable to all pipelines that start with FASTQ data
- Very high-performance
- Easy to adopt into a pipeline without being a Rust expert

## Feature Requirements
- Subsampling to an exact number of reads
- Complete handling for multiple GEM groups
- Pluggable method for assay-specific parsing of BCs and UMIs.
- "1+epsilon" passes for barcode correction: only revisit reads that need correction
- In-process aligners: BWA, STAR & interface for 'home-brew' aligners
- In-process adapter trimming (cutadapt replacement)
- Easily retain trimmed sequences: make it easy to create a fully round-trippable BAM
- State-of-the-art MARK_DUPLICATES impl
- General-purpose data-quality metrics: replace most of REPORT_BASIC & provide hooks to extend metric suite.


## Detailed Design Notes
- Design a general purpose 'read-cluster' container that stores all R1/R2/I1/I2 data & does basic
  validation of FASTQ correctness
- 
