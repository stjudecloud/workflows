# Change Log

All notable changes to this project will be documented in this file.
 
The format is based on [Keep a Changelog](http://keepachangelog.com/).

## 2025 July

### Added

- Added `util.check_fastq_and_rg_concordance` task for alignment pre-processing [#235](https://github.com/stjudecloud/workflows/pull/235).
- Added `fastp` task for FASTQ QC and read trimming [#244](https://github.com/stjudecloud/workflows/pull/244).

### Changed

- `bwa` alignment tasks now require the `read_group` argument [#244](https://github.com/stjudecloud/workflows/pull/244).
- "validation tasks" no longer have output sections [#240](https://github.com/stjudecloud/workflows/pull/240).
- "lightweight tasks" now rely on the WDL spec's default `memory` and `disks` values instead of manaully specifying something arbitrary [#240](https://github.com/stjudecloud/workflows/pull/240).
- `star.alignment` now specifies `read_one_fastqs_gz` and `read_groups` must be non-empty [#235](https://github.com/stjudecloud/workflows/pull/235).
    - immediately reverted in [#240](https://github.com/stjudecloud/workflows/pull/240).
- `star.alignment` input `read_two_fastqs_gz` is now an `Array[File]+?`, allowing it to be omitted entirely [#235](https://github.com/stjudecloud/workflows/pull/235).
    - immediately reverted in [#240](https://github.com/stjudecloud/workflows/pull/240).

### Removed

- removed interleaved FASTQ option from `samtools.bam_to_fastq` [#244](https://github.com/stjudecloud/workflows/pull/244).
- removed all uses of non-empty (`+`) array qualifier [#240](https://github.com/stjudecloud/workflows/pull/240).
- `star.alignment` no longer sorts inputs prior to passing to STAR [#235](https://github.com/stjudecloud/workflows/pull/235).
    - This means you have to ensure the FASTQs and RG records are concordant in your inputs!
- `util.get_read_groups` has been removed [#235](https://github.com/stjudecloud/workflows/pull/235).
    - see `data_structures/read_group.wdl` for an alternative.

## 2025 May

### Fixed

- Replace an object to struct coercion in `samtools.wdl` with the proper struct type [#227](https://github.com/stjudecloud/workflows/pull/227).

## 2025 February

## Added

- Added a default `prefix` calculation for `star.alignment` [#220](https://github.com/stjudecloud/workflows/pull/220).

### Changed

- `ngsderive.encoding` removed the `String inferred_encoding` output [#216](https://github.com/stjudecloud/workflows/pull/216).
- Improved the REGEX used to calculate a prefix for FASTQ input files in various tools [#220](https://github.com/stjudecloud/workflows/pull/220).
 
## 2025 January

### Changed

- `util.get_read_groups` had the param `format_for_star` reworked to a more generic `clean` parameter [#205](https://github.com/stjudecloud/workflows/pull/205).
- `star.allignment` now takes an `Array[String]` for `read_groups` instead of a single `String` [#205](https://github.com/stjudecloud/workflows/pull/205).
