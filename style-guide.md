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
  - 3 high level input sections
    - first are required inputs
    - second are optional inputs
    - last are inputs with a default value
  - within each section, inputs are ordered by variable type
    - `File`
    - `Array[*]+`
    - `Array[*]`
    - `struct`
    - `Object`
    - `Map[*, *]`
    - `Pair[*, *]`
    - `String`
    - `Boolean`
    - `Float`
    - `Int`
  - for ordering of the same compound type, drop the outermost type and recursively apply above sorting
  - once the above ordering is satisfied, it is up to the developer for final order of inputs of the same type.
  - Disallowed input names: `/^[iI]n[A-Z_]/`, `/^input/i`
- `command` blocks should be wrapped with arrows (`<<< >>>`) instead of brackets (`{ }`)
  - Certain Bash constructions cause problems with the bracket notation
- output variable names should be short but descriptive
  - Disallowed output names: `/^[oO]ut[A-Z_]/`, `/^output/i`, `/^..?$/`
- All tasks should run in a Docker container
  - linter will have an option to disable this error
- no whitespace on empty lines
- no whitespace at the end of lines
- indentation should be 4 spaces and never tab literals
- At most one empty line in a row
- End the file with a newline
- Comments on the same line as code should have 2 spaces before the `#` and one space before the comment
- Files should not mix `\n` and `\r\n` line breaks
- WDL lines should be less than 90 characters wide whenever possible
  - Exceptions would be long strings that WDL doesn't allow to be broken up
  - This restriction applies to embedded code in the `command` block as well.
