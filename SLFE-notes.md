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
- Clear separation of 'heavy-duty Rust infrastructure' and 'bioinformatics' in code.
  Ideally each 'bioinformatics' step is captured in a trait that can be implemented
  in a few simple function that don't need to know anything about advanced Rust
  features or special setups that improve performance (e.g. threading, channels,
  custom iterators, or lifetimes).

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


## Read Container 

From the FASTQs: 
- 4 read parts w/ (header, seq, qual)
- which read is which (i1 / i2 / r1 / r2)

Where the read came from:
- Gem Group (unique int)
- SAM Read Group: (flowcell, lane, sample, library)


- reads from the same Gem Group can be in different read groups
- reads from different libraries can have the sample Gem group
- any relationship constraints on sample <--> Gem group relation?

Barcode:

enum {
  Invalid(gem_group, seq)
  Valid(gem_group, seq)
}

- single sort order for all reads with or without a valid barcode.
- Sort by Invalid/Valid, Gem_group, Seq. Provide methods to get ranges
  with particular partitioning. (Tie in to shardio?)
- for shardio compat it needs to be self-contained (ie no references).
- can it be a fixed-size type (no vector allocation?)?? - otherwise it needs
  to be generic. Should be OK to fix at 16bp for the forseeable future.

- how to handle double-barcode Maverick scheme? (Ignore for now)
