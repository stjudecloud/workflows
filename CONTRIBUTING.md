# Contributing

When contributing to this repository, please first discuss the change you wish to make via issue,
email, or any other method with the owners of this repository before making a change. 

Please note we have a code of conduct, please follow it in all your interactions with the project.

## How can I start contributing?

### I don't want to write code, can I still contribute?

Sure!

We welcome bug reports. If you discover a flaw in our codebase, please review the list of open issues to ensure that it is not a duplicate. Please also attempt to debug the issue locally and ensure that it is not a configuration issue. Once you have done both, please file a new issue providing the relevant information in the issue. Please provide the exact steps to reproduce the problem, specific example(s) that demonstrate the steps, and the behavior you observe as well as the behavior you expected to observe. A copy and paste of the command and error log is always helpful (please use markdown formatting appropriately).

We also appreciate feedback on our documentation. Feel free to look over any of our `*.md` files and note any issues you find. You can also find a rendered version of our documentation at https://stjudecloud.github.io/workflows/. This is built from the embedded documentation in the various WDL files.

The maintainers reserve the right to close issues and discussions as deemed necessary as well as to delete comments and interactions within the repository.

### Your first code contribution

We encourage you to reach out to the core team prior to writing up a pull request. **This is to ensure there isn't any wasted effort or duplication of work caused by miscommunication. Failure to do so may result in the rejection of the pull request.** You can get in touch with us via the [issues][issues] or hop over to the [discussions][discussions]. We are also active on the [openwdl Slack workspace](https://openwdl.slack.com).

We encourage contributors to comment on open issues that they intend to work on to help avoid duplication of effort. If multiple individuals are interested in solving the same issue, we recommend reaching out to one another to gauge if there is potential for a collaboration.

That being said, we will not assign issues to external contributors, and commenting on an issue does not guarantee exclusive rights to work on that issue. If multiple PRs are received for the same issue, the PR that (a) most thoroughly addresses the problem being solved and (b) has the best implementation by judgement of the St. Jude Rust Labs team will be accepted in favor of the other submitted PRs.

### Review Policy

Our pull request template has an extensive checklist that must be completed prior to review. Our policy is that any PRs submitted with an incomplete checklist will not be reviewed. Part of this checklist includes ensuring that our CI checks pass. Additional guidance for satisfying the CI checks can be [found below](#the-ci-has-turned-red-how-do-i-make-it-green-again-ci-green).

Note that the maintainers reserve the right to close any submission without review for any reason.

## Expectations for WDL contributions

We have some opinionated rules and guidelines we use while writing WDL for this repository. These include:

- See `template/common-parameter-meta.txt` for common description strings.
  - If applicable, use the same parameter name and help text as the underlying tool called by the task.
- All requirement values are overridable at runtime. However, tasks should have easily configurable memory and disk space allocations.
  - See the various tasks in the template directory for possible ways to allocate resources.
    - Contributors can mix and match the available templates, copy and pasting subsections as appropriate.
    - A task may contain both statically and dynamically allocated resources.
- Multi-core tasks should *always* follow the conventions laid out in the `use_all_cores_task` example (see `template/task-examples.wdl`).
  - This is catering to cloud users, who may be allocated a machine with more cores than are specified by the `ncpu` parameter.
  - Note that future versions of WDL will likely cause a change to this convention.
    - We plan to deprecate the `ncpu` param in favor of accessing the runtime section directly (`n_cores=~{task.runtime.cpu}`).
- Output file names should *always* be determined with either the `outfile_name` parameter or the `prefix` parameter.
  - `outfile_name` should be preferred if no downstream tasks/tools rely on the file name/extension.
  - Tasks with multiple outputs should always use the `prefix` convention.
- After the input sorting rules in `sprocket lint` have been applied, follow the below rules for further sorting.
  - "sample" files come before "reference" files.
  - If present, `use_all_cores` should be the last `Boolean` in its block.
  - The `ncpu` parameter comes before inputs that allocate memory, which come before inputs that allocate disk space.
    - This block of 2-3 inputs should come after all other inputs.
- If a task uses multiple cores or is multithreaded, then at least 2 cpu should be specified.
- Use the `as` keyword sparingly; only in the case of increased readability or to avoid name collisions.
  - Prefer using `as` in the import block rather than at the task/workflow call level.
  - When using `as` to rename an invalid URI, attempt to make as few changes to the filename as possible (i.e. try not to abbreviate).
  - To disambiguate a task or workflow file from it's contents, you can respectively add the `_tasks` or `_wf` suffix in the import section.
- Whenever possible, prefer a Docker image maintained by an external source (such as BioContainers) rather than creating your own image.
- When adding a `Dockerfile` to this repository, follow the below conventions:
  - Create a directory under the `docker/` directory and choose an appropriate name (likely shared with the underlying tool). The `Dockerfile` should be nested under this new directory. Then create a `package.json` alongside the `Dockerfile`. The `package.json` file is required to contain two JSON fields (`name` and `version`). It can optionally contain a `revision` field.
  - Docker images should be versioned according to the following convention
    - The `version` should be shared with whatever underlying tool is being used
      - If no specific tool is named (e.g. the `util` image), default to SemVer. Ignore the next 3 bullet points.
    - The revision should start with zero (`0`).
      - If the Docker image gets updated, *without* updating the base tool's version, increment the number by one.
      - If the Docker image gets updated, *including* updating the base tool's version, revert back to zero.
- Any tasks which are deprecated should have a `deprecated: true` key in their `meta` section.
  - Never include a `deprecated: false` key in any production tasks. All tasks are assumed to not be deprecated unless otherwise noted.
  - In addition, there should be a `warning` key which starts with the text `**[DEPRECATED]**`.
    - No other text or explanation is required after the above text, but it can be added for further context.
  - These two conventions allow for a task's deprecated status to be communicated in multiple ways, ensuring no user misses the notice.
  - Deprecated tasks should be placed at the end of their file.
- While WDL allows embedded scripts in the `command` block sections, this repository requires scripts (e.g. R, Python) to be separate and placed in the `scripts` folder. The relevant Docker image build for your task should then include the script during the build so the task can access it. This separation of concerns improves the developer experience by improving syntax highlighting in the WDL document and enabling linting and formatting checks for the scripting languages.

## FAQs

### Can I use Artificial Intelligence (AI)?

We have found that AI, while helpful in some contexts, causes more confusion and work for all parties involved when interacting with a large, complex codebase such as the `workflows` WDL repository. To that end, no PRs including AI-generated content—whether that be generated code, generated documentation, generated discussion via GitHub comments, or any other AI generated content—will be accepted from external contributors. Any submissions deemed to be AI-generated from external contributors will be closed without review.

### What IDE should I use?

Most of this team uses VSCode with the `sprocket` extension but that preference is not hardcoded anywhere. Feel free to use any IDE you want!

### What's a good first issue?

We will try to keep a handful of [issues][issues] marked `good first issue` open and ready for new contributors.

### The CI has turned red. How do I make it green again?

The Sprocket Lint check is the most common failure. We encourage contributors to run Sprocket locally via the VSCode extension. Then address any WDL style issues before committing.

### Container Vulnerabilities

We use Snyk to scan our container images for vulnerabilities. These vulnerability lists can be quite long and therefore we prioritize Critical and High severity vulnerabilities. In addition, lower severity vulnerabilities should also be addressed when feasible. Ultimately, the maintainers of this repository are the arbiters of which vulnerabilities must be addressed to merge a pull request.

## Code of Conduct

### Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, gender identity and expression, level of experience,
nationality, personal appearance, race, religion, or sexual identity and
orientation.

### Our Standards

Examples of behavior that contributes to creating a positive environment
include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or
advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic
  address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

### Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

### Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

### Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at [INSERT EMAIL ADDRESS]. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

### Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at [http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/