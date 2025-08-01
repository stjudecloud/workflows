<p align="center">
  <a href="https://github.com/stjudecloud/workflows"><img src="./assets/workflows-banner-flowchart.jpg" width="800" title="St. Jude Cloud Workflows"></a>
  <br />
  <a href="https://actions-badge.atrox.dev/stjudecloud/workflows/goto?ref=main">
    <img alt="Build Status" src="https://img.shields.io/endpoint.svg?url=https%3A%2F%2Factions-badge.atrox.dev%2Fstjudecloud%2Fworkflows%2Fbadge%3Fref%3Dmain&style=flat" />
  </a>
  <a href="https://github.com/stjudecloud/workflows/blob/master/LICENSE.md" target="_blank">
    <img alt="License: MIT" src="https://img.shields.io/badge/License-MIT-yellow.svg" />
  </a>
</p>

> This repository contains all bioinformatics workflows used on the St. Jude Cloud project. Officially, the repository is in beta — the project is adding workflows as they are developed and put into production.

Resources requirements have been optimized to minimize failures in our computing environment, but they may not reflect the best settings for your use case. Please ensure that you tailor these parameters to fit your needs.

### 🏠 [Homepage](https://stjude.cloud)

Our [documentation](https://stjudecloud.github.io/workflows/) is generated with [Sprocket](https://sprocket.bio).

## Repository Structure

The repository is laid out as follows:

* `workflows/` - Directory containing all end-to-end bioinformatics workflows.
* `tools/` - All tools we have wrapped as individual WDL tasks.
* `data_structures/` - WDL `struct` definitions and tasks or workflows related to their construction, parsing, or validation.
* `docker/` - Dockerfiles used in our workflows. All docker images are published to the [GitHub Container Registy](https://github.com/orgs/stjudecloud/packages?repo_name=workflows) as a part of our CI and are versioned.
* `scripts/` - This directory is passed to the docker build step as a context, and contains scripts that can be copied to our docker containers.
* `tests/` - Home to all of our testing infrastructure. We use [pytest-workflow](https://pytest-workflow.readthedocs.io/en/stable/) for validating our code.
* `developer_scripts` - Home to any scripts that ease the development process.
* `bin/` - **no longer in use** Scripts used by Cromwell configuration settings. Add this to `$PATH` prior to using configurations in `conf` with Cromwell.
* `conf/` - **no longer in use** Cromwell configuration files created for various environments that we use across our team. Feel free to use/fork/suggest improvements.

## Expected FASTQ file name conventions

The tasks and workflows in this repository which have one or more FASTQ files as an input will also have a `prefix` input which will determine the filenames for any output files. The `prefix` input can be specified manually, or it can be left at the default value. The default value will attempt to strip common file suffixes from one of the input FASTQs and determine an appropriate basename to be used by all output files. That evaluation is most commonly performed using the POSIX ERE Regular Expression `(([_.][rR](?:ead)?[12])((?:[_.-][^_.-]*?)*?))?\.(fastq|fq)(\.gz)?$`. In plain English, this REGEX will at a minimum search for and remove the file extensions `.fastq` and `.fq` with or without a `.gz` GZIP extension. Additionally, if the FASTQ filename contains a "read number" signifier (`R1`/`R2`/`r1`/`r2`/`read1`/`read2`) somewhere before the FASTQ extension, that will be truncated off the basename. This means that _everything_ after the read indicator will be removed. If there is important information encoded in your filenames _between the read number and the final extension,_ we recommend you manually specify an appropriate `prefix` value.

### Examples

Every filename in the following list will have the evaluated `prefix` "`sample`" if no override value is provided.

```
sample_R1_100000.fastq.gz
sample_R2.fq
sample.fq.gz
sample.r1_100000.trimmed.fastq.gz
sample.R2_100000.trimmed-kebab.fastq.gz
sample_r1.100000.trimmed-kebab.foobar.fastq.gz
sample.Read2_100000.trimmed-kebab.fastq
sample_read1.100000.trimmed-kebab.fq
```

A FASTQ with the filename `sample.100000-kebab.foobar.fastq.gz` would have a default `prefix` value of "`sample.100000-kebab.foobar`".

The following filenames will not result in _any_ trimming of the filename, and should likely be either renamed or have a manually specified `prefix`:

```
sample_R_one.FASTQ.gz
sampleR1.Fq
sample_read_two.fq.zip
```

## Bootstrap guide

This repository implements workflows using the Workflow Description Language (WDL). If unfamiliar with WDL, a short overview is available in the [WDL spec](https://github.com/openwdl/wdl/blob/main/versions/1.1/SPEC.md#introduction).

The workflows and tasks in this repository should require minimal set-up and configuration before you're ready to run. You don't even need to clone the repo! The bare minimum requirements are a locally installed WDL runner and an internet connection.

The exact steps for installation, configuration, and execution are going to depend on you environment and preferred engine. There are a variety of WDL engines you could use, though our team prefers [miniwdl](https://github.com/chanzuckerberg/miniwdl). We also make use of the [`miniwdl-lsf` plugin](https://pypi.org/project/miniwdl-lsf/) for running on our LSF cluster.

Most WDL runners are capable of running a WDL file from a URL. This is how we most commonly execute our workflows and tasks. The below command is a mock example of of what could be used to submit a run of our rnaseq-standard workflow using `miniwdl`:

```bash
miniwdl run --verbose --input inputs.json https://raw.githubusercontent.com/stjudecloud/workflows/rnaseq-standard/v3.0.1/workflows/rnaseq/rnaseq-standard.wdl
```

For an introduction to WDL, there are many guides, one of which is [from Terra](https://support.terra.bio/hc/en-us/articles/360037117492-Overview-Getting-started-with-WDL).

## Author

👤 **St. Jude Cloud Team**

* Website: https://stjude.cloud
* Github: [@stjudecloud](https://github.com/stjudecloud)
* Twitter: [@StJudeResearch](https://twitter.com/StJudeResearch)

## Tests

Every task in this repository is covered by at least one test (see all of our tests in `tests/tools/`). These are run using [pytest-workflow](https://pytest-workflow.readthedocs.io/en/stable/).

The command for running our tests should be executed at the root of the repo: `python -m pytest --kwdof --git-aware`

## 🤝 Contributing

Contributions, issues and feature requests are welcome!<br />Feel free to check the [issues page](https://github.com/stjudecloud/workflows/issues). You can also take a look at the [contributing guide](https://github.com/stjudecloud/workflows/blob/master/CONTRIBUTING.md).

## Links worth checking out

[The OpenWDL GitHub](https://github.com/openwdl)

Our preferred WDL runner: [miniwdl](https://github.com/chanzuckerberg/miniwdl)

Most of our tasks are run inside a [BioContainers image](https://github.com/BioContainers/containers)

Our tasks are validated using [pytest-workflow](https://pytest-workflow.readthedocs.io/en/stable/)

## 📝 License

Copyright © 2020-Present [St. Jude Cloud Team](https://github.com/stjudecloud).<br />
This project is [MIT](https://github.com/stjudecloud/workflows/blob/master/LICENSE.md) licensed.
