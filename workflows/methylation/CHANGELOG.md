# Change Log

All notable changes to this project will be documented in this file.
 
The format is based on [Keep a Changelog](http://keepachangelog.com/).
 
## 2025 September

### Added

- Now computes a per-sample [detection value](https://www.rdocumentation.org/packages/minfi/versions/1.18.4/topics/detectionP) for each probe and provides functionality to filter based on the fraction of a cohort with p-values exceeding a threshold [#262](https://github.com/stjudecloud/workflows/pull/262)

## 2025 August

### Changed

- Beta value matrix now retains sample names as column names [#259](https://github.com/stjudecloud/workflows/pull/259)

## 2025 July

### Changed

- "lightweight tasks" now rely on the WDL spec's default `memory` and `disks` values instead of manaully specifying something arbitrary [#240](https://github.com/stjudecloud/workflows/pull/240).

### Removed

- removed all uses of non-empty (`+`) array qualifier [#240](https://github.com/stjudecloud/workflows/pull/240).
 
## 2025 May

### Changed

 * Allow nested inputs [#229](https://github.com/stjudecloud/workflows/pull/229)
