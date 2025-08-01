# Change Log

All notable changes to this project will be documented in this file.
 
The format is based on [Keep a Changelog](http://keepachangelog.com/).
 
## 2025 July

### Added

- Added `fastp` QC and optional read trimming (also from `fastp`) [#247](https://github.com/stjudecloud/workflows/pull/247).
- Added read group validation [#235](https://github.com/stjudecloud/workflows/pull/235).

### Fixed

- FASTQ entrypoint now specifies that FASTQ and read group record inputs must be non-empty [#235](https://github.com/stjudecloud/workflows/pull/235).

### Changed

- "validation tasks" no longer have output sections [#240](https://github.com/stjudecloud/workflows/pull/240).
- "lightweight tasks" now rely on the WDL spec's default `memory` and `disks` values instead of manaully specifying something arbitrary [#240](https://github.com/stjudecloud/workflows/pull/240).
- Added `after` clauses to many calls to prevent wasted compute on validation failures [#235](https://github.com/stjudecloud/workflows/pull/235).

### Removed

- removed all uses of non-empty (`+`) array qualifier [#240](https://github.com/stjudecloud/workflows/pull/240).

## 2025 February

## Added

- Added a default `prefix` calculation for `dnaseq-standard-fastq` [#220](https://github.com/stjudecloud/workflows/pull/220).