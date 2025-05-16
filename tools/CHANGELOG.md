# Change Log

All notable changes to this project will be documented in this file.
 
The format is based on [Keep a Changelog](http://keepachangelog.com/).

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
