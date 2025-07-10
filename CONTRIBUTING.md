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

## FAQs

### Can I use Artificial Intelligence (AI)?

We have found that AI, while helpful in some contexts, causes more confusion and work for all parties involved when interacting with a large, complex codebase such as the `workflows` WDL repository. To that end, no PRs including AI-generated content—whether that be generated code, generated documentation, generated discussion via GitHub comments, or any other AI generated content—will be accepted from external contributors. Any submissions deemed to be AI-generated from external contributors will be closed without review.

### What IDE should I use?

Most of this team uses VSCode with the `sprocket` extension but that preference is not hardcoded anywhere. Feel free to use any IDE you want!

### What's a good first issue?

We will try to keep a handful of [issues][issues] marked `good first issue` open and ready for new contributors.

### The CI has turned red. How do I make it green again?

The Sprocket Lint check is the most common failure. We encourage contributors to run Sprocket locally via the VSCode extension. Then address any WDL style issues before committing.

The Docker Pull Check also fails when an image version is incremented, but the corresponding change has not made it to the `main` branch to build the docker image. These failures can be ignored if the pull request is responsible for creating the new docker image version.

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