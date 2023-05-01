# WDL Style Guide

All rules below should be followed by contributors to this repo. Contributors should also follow the less specific rules outlined in `style-guide.md` Pull Requests which do not conform to these specifications will be asked to change.

## Rules

- All WDL should be written in v1.0
- All inputs must have a corresponding `parameter_meta` entry
  - These texts should be copy and pasted from other tasks with the same input when possible
    - see `template/common-parameter-meta.txt` for common description strings.
  - If applicable, use the same parameter name, help string, and parameter ordering as the underlying tool called by the task
- TODO ordering of inputs more specifically than `style-guide.md`
  - e.g. RAM allocations before disk space allocations
- Tasks with string parameters for which a limited number of choices are valid, must be documented following the template in `string_choices_task` (see `template/task-templates.wdl`)
  - they should also fail quickly with an informative error message if an invalid input is provided
    - In most cases, just passing the parameter to the underlying tool should produce a satisfactory error, but this must be checked for each task
- All tasks must have configurable memory and disk space allocations
  - see the various tasks in the template directory for possible ways to allocate resources
    - Contributors can mix and match the available templates, copy and pasting subsections as appropriate
    - It is allowed to have one resource allocated dynamically, and another allocated statically in the same task.
    - It is *not* allowed to have a resource which can be allocated *either* statically or dynamically.
      - This is technically feasible, but is too complicated for maintenance and end-users.
      - e.g. `memory_gb` and `modify_memory_gb` cannot be present in the same task.
- All tasks and workflows should have a `max_retries` input.
  - This should be defaulted to `1` for nearly all tasks
  - Some tasks are particularly error prone and can have a higher default `max_retries`
  - **rule specific to workflows:** `max_retries` should be an optional `Int?`
    - This allows each task to have it's own specific default `max_retries`
      - If a user does not supply `max_retries`, those task level defaults will get used
      - If a user does supply `max_retries`, it should override the default for *every* task called
- multi-core tasks should *always* follow the conventions laid out in the `detect_nproc_task` example (see `template/task-templates.wdl`)
  - this is catering to cloud users, who may be allocated a machine with more cores than are specified by the `ncpu` parameter
- output file names should *always* be determined with either the `outfile_name` parameter or the `prefix` parameter.
  - `outfile_name` should be preferred if no downstream tasks/tools rely on the file name/extension
  - tasks with multiple outputs should always use the `prefix` convention
- All tasks should run in a Docker container
  - whenever possible, prefer an image maintained by an external source (such as BioContainers) rather than creating your own image
  - general purpose tasks can use the `util` image maintained in this repo
