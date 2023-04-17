# WDL Style Guide

All rules below should be followed by contributers to this repo. Pull Requests which do not conform to these specifications will be asked to change.

## Rules

- All WDL should be written in v1.0
- The following sections must be present and in this order: `meta`, `parameter_meta`, `input`, `command`, `output`, `runtime`
- The `meta` section should have a `description` of the task
  - `meta` should *not* include `author` or `email` keys
- All inputs must have a corresponding `parameter_meta` entry
  - These texts should be copy and pasted from other tasks with the same input when possible
    - see `template/common-parameter-meta.txt` for common description strings.
  - Inputs and parameter meta must be in the same order
  - That order should be logical
    - No hard and fast rules, but look to similar tasks/workflows for an example of proper ordering
- Tasks with string parameters for which a limited number of choices are valid, must be documented following the template in `string_choices_task` (see `template/task-templates.wdl`)
  - they should also fail quickly with an informative error message if an invalid input is provided
    - In most cases, just passing the parameter to the underlying tool should produce a satisfactory error, but this must be checked for each task
- All tasks must have configurable memory and disk space allocations
  - see the various tasks in the template directory for possible ways to allocate resources
    - Contributors can mix and match the available templates, copy and pasting subsections as appropriate
- All tasks should have a `max_retries` input.
  - This should be defaulted to `1` for nearly all tasks
  - Some tasks are particularly error prone and can have a higher default `max_retries`
  - **rule specific to workflows:** `max_retries` should be an optional `Int?`
    - This allows each task to have it's own specific default `max_retries`
      - If a user does not supply `max_retries`, those task level defaults will get used
      - If a user does supply `max_retries`, it should override the default for *every* task called
- multi-core tasks should always follow the conventions laid out in the `detect_nproc_task` example (see `template/task-templates.wdl`)
- `command` blocks should be wrapped with arrows (`<<< >>>`) instead of brackets (`{ }`)
  - Certain Bash constructions cause problems with the bracket notation
- output file names should *always* be determined with either the `outfile_name` parameter or the `prefix` parameter.
  - TODO: address tasks with multiple outputs
- output variable names should be short but descriptive
- All tasks should run in a Docker container
- no whitespace on empty lines
- WDL lines should be less than 100 characters wide whenever possible
  - Exceptions would be long strings that WDL doesn't allow to be broken up
  - TODO: Do we enforce this within the `command` block? Or do we allow Bash to exceed 100 characters?
