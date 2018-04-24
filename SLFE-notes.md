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


## Components & current status
- Rust metrics framework - pretty complete - (where is this?)
- Shardio - v2 need API review, switch existing clients over, polish - (https://github.com/10XGenomics/rust-shardio/tree/pmarks/bam-gen)
- Read / GemGroup / ReadGroup / Barcode / UMI trait scheme - need lots more design work - (fastq/src/new./rs)
- cutadapt clone - needs lots of work to cover cases used by Agora, GEX, VDJ - (see fastq/src/cutadapt.rs) 
- rust-bwa - in-process BWA wrapper. in a working state, needs stress-testing : https://github.com/10XGenomics/rust-bwa/
- rust-STAR - in-process STAR wrapper. not really started.
- Better ergonomics for writing Rust stages? (see https://github.com/10XDev/TDA_rust/blob/master/src/cmd_sort_fastq.rs#L613)
  - Can we do type-driven interface & get rid of JSON handling?
  - Martian code-gen?
  - What would the simplest possible stage definition look like?
