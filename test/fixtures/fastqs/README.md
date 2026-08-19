# Input Files used for testing

All the files in this directory were either randomly generated or sourced from publicly available genomic data. Many of the files have been downsampled to make them a manageable size for storing in a Git repository (using `git-lfs`) and to enable quick test execution in our CI. For the most part, these files are not representative of any meaningful scientific data and should be considered as mocked data only suitable for the purpose of generic testing.

The following list is sorted alphabetically:

## random10k.r1.fq.gz

10,000 randomly generated read1 FASTQ reads. Generated with the [`fq` tool](https://github.com/stjude-rust-labs/fq). Gzipped.

## random10k.r2.fq.gz

10,000 randomly generated read2 FASTQ reads. Generated with the [`fq` tool](https://github.com/stjude-rust-labs/fq). Gzipped.

## test_R1.fq.gz

10,000 reads in FASTQ format. Can be used by itself to represent a Single-End sample or used with `test_R2.fq.gz` to represent a Paired-End sample. Generated with `fq generate`.

## test_R2.fq.gz

10,000 reads in FASTQ format. Can be used by itself to represent a Single-End sample or used with `test_R1.fq.gz` to represent a Paired-End sample. Generated with `fq generate`.

## test_variants_R1.fq.gz

10,000 reads in FASTQ format. Can be used by itself to represent a Single-End sample or used with `test_variants.R2.fq.gz` to represent a Paired-End sample. Generated with `fq generate` on `chrY` and `chrM`. Simulated from a reference sequence containing variants.

## test_variants_R2.fq.gz

10,000 reads in FASTQ format. Can be used by itself to represent a Single-End sample or used with `test_variants.R1.fq.gz` to represent a Paired-End sample. Generated with `fq generate` on `chrY` and `chrM`. Simulated from a reference sequence containing variants.