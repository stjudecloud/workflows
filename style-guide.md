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
- All import statements should follow the WDL version declaration (with one empty line between the version and the first import statement)
- Import statements should be sorted by the lexicographical ordering of each entire line
  - No extra white space allowed between symbols or lines
- For workflows, the following sections must be present and in this order: `meta`, `parameter_meta`, `input`, `output`
  - `input`, `parameter_meta`, and `output` are technically optional in the WDL spec, though it is discouraged to write workflows in this manner.
    - The linter will not enforce the presence of these sections
      - If `input` is present, the linter will enforce the presence of `parameter_meta`
- For tasks, the following sections must be present and in this order: `meta`, `parameter_meta`, `input`, `command`, `output`, `runtime`
  - `input`, `parameter_meta`, and `output` are technically optional in the WDL spec, though it is discouraged to write tasks in this manner.
    - The linter will not enforce the presence of these sections
      - If `input` is present, the linter will enforce the presence of `parameter_meta`
- The `meta` section should have a `description` of the task or workflow
  - The `description` should be in active voice, beginning the first sentence with a verb
    - Each task/workflow is _doing_ something. The first sentence should be a succinct description of what that "something" is.
- The `meta` section should have an `outputs` key and keys with descriptions for each output of the task/workflow
- Additional `meta` entries are allowed (such as `author` or `email` keys)
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
- There should be an empty line between workflow/task sections and any intermediate declarations
- End the file with a newline
- Comments on the same line as code should have 2 spaces before the `#` and one space before the comment
- Files should not mix `\n` and `\r\n` line breaks
- WDL lines should be less than 90 characters wide whenever possible
  - Exceptions would be long strings that WDL doesn't allow to be broken up
  - This restriction applies to embedded code in the `command` block as well.
