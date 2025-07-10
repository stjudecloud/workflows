# Change Log

All notable changes to this project will be documented in this file.
 
The format is based on [Keep a Changelog](http://keepachangelog.com/).
 
## 2025 July

### Changed

- `read_group_to_string` task now named `inner_read_group_to_string`. This task should not be called by end users [#235](https://github.com/stjudecloud/workflows/pull/235).

### Added

- `read_group_to_string` workflow, which validates `ReadGroup` structs prior to converting them to strings [#235](https://github.com/stjudecloud/workflows/pull/235).
    - This should be used instead of the `inner_read_group_to_string` task.
