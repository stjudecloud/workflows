# WDL Style Guide

All rules below should be followed by contributors to this repo. Pull Requests which do not conform to these specifications will be asked to change. This small set of rules is in addition to the lint rules enforced by Sprocket. 

## Rules

- The `description` should be in active voice, beginning the first sentence with a verb
  - Each task/workflow is _doing_ something. The first sentence should be a succinct description of what that "something" is.
- If documenting a workflow, task, input, or output and you need to be more verbose than is appropriate in a `description:` field, you may include _in addition_ an `external_help:` key with a URL
  - `external_help` is _not_ a substitute for internal documentation, although it may allow the internal documentation to be briefer
- All tasks should run in a persistently versioned container
  - This ensures reproducibility across time and environments
- Any tasks which are deprecated should have a `deprecated: true` key in their `meta` section
  - It is allowed (but redundant and discouraged) to include a `deprecated: false` key in any production tasks. All tasks are assumed to not be deprecated unless otherwise noted.
  - In addition, the `description` key of deprecated tasks should start with `**[DEPRECATED]**`
    - These two rules allow for a task's deprecated status to be communicated in multiple ways, ensuring no user misses the notice
- Deprecated tasks should be placed at the end of their file
