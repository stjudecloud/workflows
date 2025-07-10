# Change Log

All notable changes to this project will be documented in this file.
 
The format is based on [Keep a Changelog](http://keepachangelog.com/).
 
## NOTICE

Any change to the `scripts/` directory should be accompanied by version increases in the `docker/` directory! If you are editing this file, please ensure these changes propagate!

## 2025 July

### Added

- Added `util/check_FQs_and_RGs.py` script [#235](https://github.com/stjudecloud/workflows/pull/235).
    - This script is a refactor of `star/sort_star_input.py` which is generalized for cases other than STAR.

### Removed

- Removed `star/sort_star_input.py` [#235](https://github.com/stjudecloud/workflows/pull/235).
    - As a consequence, the custom STAR docker image is no more.

## 2025 February

### Added

- Added `calc_global_phred_scores.py` and `calc_gene_lengths.py` scripts [#216](https://github.com/stjudecloud/workflows/pull/216).

## 2025 January

### Added

- Added `black` formatting and `pyright` validation for Python scripts [#201](https://github.com/stjudecloud/workflows/pull/201).

### Fixed

- Fixed processing of STAR inputs [#202](https://github.com/stjudecloud/workflows/pull/202).
- _Actually_ fixed processing of STAR inputs [#205](https://github.com/stjudecloud/workflows/pull/205).
