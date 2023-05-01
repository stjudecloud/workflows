# WDL Style Guide

All rules below should be followed by contributors to this repo. Pull Requests which do not conform to these specifications will be asked to change.

These rules might also be enforced by a yet-to-be-written linter.

## Rules

- The following sections must be present and in this order: `meta`, `parameter_meta`, `input`, `command`, `output`, `runtime`
- The `meta` section should have a `description` of the task
  - `meta` should *not* include `author` or `email` keys
- All inputs must have a corresponding `parameter_meta` entry
  - These texts should be copy and pasted from other tasks with the same input when possible
  - Inputs and parameter meta must be in the same order
  - TODO: rework this section
  - That order should be logical
    - look to similar tasks/workflows for an example of proper ordering
    - Generally speaking:
      - Required inputs at the top
      - "sample" files before reference files
      - Resource configuration at the bottom
        - memory allocations come before disk space allocations
      - `max_retries` should be last
      - If applicable, `detect_nproc` immediately before `max_retries`
- Tasks with string parameters for which a limited number of choices are valid, must be documented following the template in `string_choices_task` (see `template/task-templates.wdl`)
  - TODO: this can't be enforced by a linter. Do we keep it in the guide?
  - they should also fail quickly with an informative error message if an invalid input is provided
    - In most cases, just passing the parameter to the underlying tool should produce a satisfactory error, but this must be checked for each task
- All tasks must have configurable memory and disk space allocations
  - TODO: only enforceable if we require specific var naming. Do we want to enforce specific variable names? Keep this in the guide?
  - see the various tasks in the template directory for possible ways to allocate resources
    - Contributors can mix and match the available templates, copy and pasting subsections as appropriate
    - It is allowed to have one resource allocated dynamically, and another allocated statically in the same task.
    - It is *not* allowed to have a resource which can be allocated *either* statically or dynamically.
      - This is technically feasible, but is too complicated for maintenance and end-users.
      - e.g. `memory_gb` and `modify_memory_gb` cannot be present in the same task.
- All tasks and workflows should have a `max_retries` input.
  - TODO: again, only enforeceable with specific var names. Not sure if we want to get into that business. Keep it in the guide?
  - This should be defaulted to `1` for nearly all tasks
  - Some tasks are particularly error prone and can have a higher default `max_retries`
  - **rule specific to workflows:** `max_retries` should be an optional `Int?`
    - This allows each task to have it's own specific default `max_retries`
      - If a user does not supply `max_retries`, those task level defaults will get used
      - If a user does supply `max_retries`, it should override the default for *every* task called
- multi-core tasks should *always* follow the conventions laid out in the `detect_nproc_task` example (see `template/task-templates.wdl`)
  - this is catering to cloud users, who may be allocated a machine with more cores than are specified by the `ncpu` parameter
  - TODO: this would be tough to validate with a linter. It is definitely a best practice, but IMO doesn't belong in the guide
- `command` blocks should be wrapped with arrows (`<<< >>>`) instead of brackets (`{ }`)
  - Certain Bash constructions cause problems with the bracket notation
- output file names should *always* be determined with either the `outfile_name` parameter or the `prefix` parameter.
  - `outfile_name` should be preferred if no downstream tasks/tools rely on the file name/extension
  - tasks with multiple outputs should always use the `prefix` convention
  - TODO: IMO another "best practice" that shouldn't be enforced by a linter. Remove from guide?
- output variable names should be short but descriptive
  - TODO: This I'm more inclined to enforce with the linter. I think disallowing specific output names is better than forcing users into using specific ones.
  - Disallowed output names: `out*`, any name less than 3 characters
    - TODO: is this sufficient?
- All tasks should run in a Docker container
  - TODO: is there a valid reason not to use a Docker container? Should this be enforced by the linter?
- no whitespace on empty lines
- no whitespace at the end of lines
- WDL lines should be less than 100 characters wide whenever possible
  - Exceptions would be long strings that WDL doesn't allow to be broken up
  - This restriction applies to embedded code in the `command` block as well.
    - TODO: I'm still not sure about this decision. I have reservations applying this rule to the command block. I've found a small handful of cases where I think the most readable formatting extends over 100 characters. They are rare, but exist frequently enough I'm not sure I like this rule.
