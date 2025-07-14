# Change Log

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/).

## 2025 July

### Changed

- "lightweight tasks" now rely on the WDL spec's default `memory` and `disks` values instead of manaully specifying something arbitrary [#240](https://github.com/stjudecloud/workflows/pull/240).

### Removed

- removed all uses of non-empty (`+`) array qualifier [#240](https://github.com/stjudecloud/workflows/pull/240).
 
## 2025 May

### Fixed

- Replace an object to struct coercion with the proper struct type [#229](https://github.com/stjudecloud/workflows/pull/229).
