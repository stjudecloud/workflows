# Change Log

All notable changes to this project will be documented in this file.
 
The format is based on [Keep a Changelog](http://keepachangelog.com/).

## 2025 September

### Changed

- `calc_tpm` task in `htseq.wdl` had the `gene_lengths` input renamed to `feature_lengths` [#269](https://github.com/stjudecloud/workflows/pull/269)
- Added `rm` statements to clean up temporary files at the end of `bwa` and `star` tasks [#268](https://github.com/stjudecloud/workflows/pull/268)
- Bumped `fq` versions from `0.11.0` to `0.12.0` [#261](https://github.com/stjudecloud/workflows/pull/261).
- Ported the python heredoc used in `htseq.calc_tpm` to a full script [#266](https://github.com/stjudecloud/workflows/pull/266)
- Updated to the latest MultiQC (v1.31) which stabilized the `.parquet` output format file [#264](https://github.com/stjudecloud/workflows/pull/264).

## 2025 August

### Changed

- `util.split_string` now requires an explicitly provided `delimiter` input (instead of defaulting to `" , "`) [#260](https://github.com/stjudecloud/workflows/pull/260).
- Scaled back `bwa.bwa_aln_pe` disk usage by 5gb [#260](https://github.com/stjudecloud/workflows/pull/260).
- Updated to the latest version of MultiQC (v1.30) and made some API changes in the process [#258](https://github.com/stjudecloud/workflows/pull/258).

### Fixed

- Upstream bug in Sprocket patched locally [#251](https://github.com/stjudecloud/workflows/pull/251)
    - Upstream bug [here](https://github.com/stjude-rust-labs/wdl/issues/574)
    - Local fix was moving a private declaration in `samtools.merge_sam_files` to the command body 

### Removed

- Removed `util.qc_summary` which has been out of date from `quality-check-standard` for a very long time [#251](https://github.com/stjudecloud/workflows/pull/251)

## 2025 July

### Added

- Added `util.check_fastq_and_rg_concordance` task for alignment pre-processing [#235](https://github.com/stjudecloud/workflows/pull/235).
- Added `fastp` task for FASTQ QC and read trimming [#244](https://github.com/stjudecloud/workflows/pull/244).

### Changed

- `util.check_fastq_and_rg_concordance` task no longer accepts tab literal delimiters. Use `\t` instead [#247](https://github.com/stjudecloud/workflows/pull/247).
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
