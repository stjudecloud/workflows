# WDL Style Guide

All rules below should be followed by contributors to this repo. Pull Requests which do not conform to these specifications will be asked to change.

These rules might also be enforced by a yet-to-be-written linter.

## Rules

- The first line of each file should be `## # <title of file>`
  - This ensures that if using WDLdoc for your documentation, it will have a proper Markdown title at the top
  - If missing, the linter would autogenerate a title based on the filename
    - This behavior could be disabled
- Files should have an [SPDX formatted](https://spdx.github.io/spdx-spec/v2.3/using-SPDX-short-identifiers-in-source-files/) Copyright and License information
  - The linter would autogenerate these fields for the user to fill in
    - or they could be filled in automatically if specified in the configuration/parameters
    - This behavior could be disabled
- The following sections must be present and in this order for all workflows: `meta`, `parameter_meta`, `input`, `output`
- The following sections must be present and in this order for all tasks: `meta`, `parameter_meta`, `input`, `command`, `output`, `runtime`
- The `meta` section should have a `description` of the task or workflow
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
- Files with workflows in them should document the outputs in a WDLdoc header
  - The linter will create a template to be filled in
    - This behavior could be disabled
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
