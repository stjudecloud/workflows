# WDL Best Practices

All rules below should be followed by contributors to this repo. Contributors should also follow the rules outlined in `style-guide.md`. Pull Requests which do not conform to these specifications will be asked to change.

## Rules

- All WDL should be written in v1.0
  - This is to enable the broadest base of support of your workflows
  - This recommendation is subject to change come broader support for WDL v1.1
- Variables should be in "snake_case"
- See `template/common-parameter-meta.txt` for common description strings.
  - If applicable, use the same parameter name, help string, and parameter ordering as the underlying tool called by the task
- Check all assumptions made about inputs before beginning long running executions that will fail if assumptions don't hold
  - Common examples of assumptions that should be checked: valid `String` choice, mutually exclusive parameters, missing optional file for selected parameters, filename extensions
  - This can commonly be handled by a `parse_input` task (defined in the same file as the workflow in question)
    - When possible, avoid passing in entire files to the `parse_input` task. Coerce files to `Boolean`s or `String`s to avoid unnecessary disk space usage
- Tasks with string parameters for which a limited number of choices are valid, must be documented following the template in `string_choices_task` (see `template/task-templates.wdl`)
  - they should also fail quickly with an informative error message if an invalid input is provided
    - In most cases, just passing the parameter to the underlying tool should produce a satisfactory error, but this must be checked for each task
  - While redundant, it is still best practice to validate these strings in the `parse_input` task of any workflow which calls the task
    - This ensures the workflow will fail as fast as possible to save users time and resources
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
- After the input sorting rules in `style-guide.md` have been applied, follow the below rules for further sorting.
  - "sample" files come before "reference" files
  - If present, `detect_nproc` should be the last `Boolean` in its block
  - the `ncpu` parameter comes before inputs that allocate memory, which come before inputs that allocate disk space, which come before `max_retries`
    - This block of 3-4 inputs should come after all other inputs.
- All tasks should have an output
  - This may be a hardcoded "dummy" output such as `String check = "passed"`
  - This ensures the task can be cached by runners. Tasks without outputs may be required to rerun on the same input due to a cache miss.
- Whenever possible, prefer a Docker image maintained by an external source (such as BioContainers) rather than creating your own image
- When adding a Dockerfile to this repository, follow the below conventions
  - The `Dockerfile` should be nested under the `docker/` directory, a folder with a name for the image (in most cases the name of the primary tool), and finally a folder named after the version being built.
  - Docker images should be versioned according to the following convention
    - Start with the version of whatever tool is named in the path to the `Dockerfile`
      - If no specific tool is named (e.g. the `util` image), default to SemVer. Ignore the next 3 bullet points.
    - Followed by a dash-zero (`-0`)
      - If the Docker image gets updated, *without* updating the base tool's version, increment the number after the dash (`-`) by one
      - If the Docker image gets updated, *including* updating the base tool's version, revert back to a dash-zero (`-0`)
- general purpose tasks can use the `util` image maintained in this repo