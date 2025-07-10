# Change Log

All notable changes to this project will be documented in this file.
 
The format is based on [Keep a Changelog](http://keepachangelog.com/).
 
## 2025 July

### Added

- Added read group validation [#235](https://github.com/stjudecloud/workflows/pull/235).

### Fixed

- FASTQ entrypoint now specifies that FASTQ and read group record inputs must be non-empty [#235](https://github.com/stjudecloud/workflows/pull/235).

### Changed

- Added `after` clauses to many calls to prevent wasted compute on validation failures [#235](https://github.com/stjudecloud/workflows/pull/235).

## 2025 February

## Added

- Added a default `prefix` calculation for `dnaseq-standard-fastq` [#220](https://github.com/stjudecloud/workflows/pull/220).