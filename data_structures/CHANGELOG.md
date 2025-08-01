# Change Log

All notable changes to this project will be documented in this file.
 
The format is based on [Keep a Changelog](http://keepachangelog.com/).

## 2025 July

### Changed

- `read_group_to_string` now inserts tab escape sequences (`\t`) instead of tab literals when `format_as_sam_record = true` [#247](https://github.com/stjudecloud/workflows/pull/247).
- "validation tasks" no longer have output sections [#240](https://github.com/stjudecloud/workflows/pull/240).
- "lightweight tasks" now rely on the WDL spec's default `memory` and `disks` values instead of manaully specifying something arbitrary [#240](https://github.com/stjudecloud/workflows/pull/240).
- `validate_string_is_12bit_oct_dec_or_hex` renamed to `validate_string_is_12bit_int` [#240](https://github.com/stjudecloud/workflows/pull/240).
- `read_group_to_string` task now named `inner_read_group_to_string`. This task should not be called by end users [#235](https://github.com/stjudecloud/workflows/pull/235).

### Added

- `read_group_to_string` workflow, which validates `ReadGroup` structs prior to converting them to strings [#235](https://github.com/stjudecloud/workflows/pull/235).
    - This should be used instead of the `inner_read_group_to_string` task.
