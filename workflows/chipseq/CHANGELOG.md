# Change Log

All notable changes to this project will be documented in this file.
 
The format is based on [Keep a Changelog](http://keepachangelog.com/).
 
## 2025 July

### Added

- Added `fastp` QC and optional read trimming (also from `fastp`) [#247](https://github.com/stjudecloud/workflows/pull/247).
- Added read group validation [#235](https://github.com/stjudecloud/workflows/pull/235).

### Changed

- Added `after` clauses to many calls to prevent wasted compute on validation failures [#235](https://github.com/stjudecloud/workflows/pull/235).

### Removed

- Removed `ngsderive readlen` call that wasn't being output [#235](https://github.com/stjudecloud/workflows/pull/235).
