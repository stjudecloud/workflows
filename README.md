<p align="center">
  <a href="https://github.com/stjudecloud/workflows"><img src="./docs/workflows-banner-flowchart.jpg" width="800" title="St. Jude Cloud Workflows"></a>
  <a href="https://actions-badge.atrox.dev/stjudecloud/workflows/goto">
    <img alt="Build Status" src="https://img.shields.io/endpoint.svg?url=https%3A%2F%2Factions-badge.atrox.dev%2Fstjudecloud%2Fworkflows%2Fbadge&style=flat" />
  </a>
  <a href="https://stjude.cloud/docs" target="_blank">
    <img alt="Documentation" src="https://img.shields.io/badge/documentation-yes-brightgreen.svg" />
  </a>
  <a href="https://github.com/stjudecloud/workflows/blob/master/LICENSE.md" target="_blank">
    <img alt="License: MIT" src="https://img.shields.io/badge/License-MIT-yellow.svg" />
  </a>
</p>

> This repository contains all bioinformatics workflows used on the St. Jude Cloud project. Officially, the repository is in beta — the project is adding workflows as they are developed and put into production.

### 🏠 [Homepage](https://stjude.cloud)

## Getting Started

At the time of writing, all workflows are written in [WDL][wdl] and are tested
using [Cromwell][cromwell]. We use [Oliver][oliver] to easily interact with the
Cromwell server to perform various tasks. Although we do not test outside of Cromwell, we
expect that the workflows will work just as well using other runners.

The easiest way to get started is to install [bioconda][bioconda] and the run the following commands:

```bash
conda create -n workflows-dev -c conda-forge cromwell miniwdl -y
conda activate workflows-dev
git clone git@github.com:stjudecloud/workflows.git
cd workflows
```

Any of the workflows in [the workflows](https://github.com/stjudecloud/workflows/tree/master/workflows) folder is a good place to start, e.g.

```bash
cromwell run workflows/reference/bootstrap-reference.wdl --inputs workflows/reference/inputs.json
```

## Repository Structure

The repository is laid out as follows:

* `bin` - Scripts used by Cromwell configuration settings. Add this to `$PATH` prior to using configurations  in `conf` with Cromwell.
* `conf` - Cromwell configuration files created for various environments that we use across our team. Feel free to use/fork/suggest improvements.
* `docker` - Dockerfiles used in our workflows. All docker images are published to [Docker Hub](https://hub.docker.com/u/stjudecloud) as a part of our CI and are versioned.
* `tools` - All tools we have wrapped as individual WDL tasks.
* `workflows` - Directory containing all end-to-end bioinformatics workflows.

## Workflows Available

The current workflows exist in this repo with the following statuses:

| Name                               | Version         | Description                                                                                                                                           | Specification                                                                                         | Workflow                                                                                                                                                                                                                                                                                                                                 | Status                                                                                                              |
| ---------------------------------- | --------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------- |
| Standard RNA-Seq                   | v3.0.0          | Standard RNA-Seq harmonization pipeline.                                                                                                              | [Specification](https://stjudecloud.github.io/rfcs/0001-rnaseq-workflow-v2.0.html)                    | [Realign BAM Workflow](./workflows/rnaseq/rnaseq-standard.wdl), [FastQ Workflow](./workflows/rnaseq/rnaseq-standard-fastq.wdl)                                                                                                                                                                                                           | ![In Production](https://img.shields.io/static/v1?label=Status&message=Production&color=green&style=flat-square)    |
| Build STAR References              | N/A             | Build [STAR aligner](https://github.com/alexdobin/STAR) reference files used in RNA-Seq Standard harmonization pipelines.                             | None                                                                                                  | [Workflow](./workflows/rnaseq/rnaseq-star-db-build.wdl)                                                                                                                                                                                                                                                                                  | ![In Production](https://img.shields.io/static/v1?label=Status&message=Production&color=green&style=flat-square)    |
| Quality Check Standard             | v1.0.0          | Perform ~10 different QC analyses on a BAM file and compile the results using [MultiQC](https://multiqc.info/).                                       | [Specification](https://rfcs.stjude.cloud/branches/rfcs/qc-workflow/0002-quality-check-workflow.html) | [Workflow](./workflows/qc/quality-check-standard.wdl)                                                                                                                                                                                                                                                                                    | ![In Production](https://img.shields.io/static/v1?label=Status&message=Production&color=green&style=flat-square)    |
| Build FastQ Screen References      | N/A             | Build references used in WGS/WES Quality Check pipeline for running [FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/). | None                                                                                                  | [Workflow](./workflows/qc/make-qc-reference.wdl)                                                                                                                                                                                                                                                                                         | ![In Production](https://img.shields.io/static/v1?label=Status&message=Production&color=green&style=flat-square)    |
| RNA-Seq Expresssion Classification | v1.0.0(?)       | Plots given samples on a labelled t-SNE graph built using expression data of (?) samples.                                                             | None(?)                                                                                               | [Realign BAM Worfklow](./workflows/rnaseq_expression_classification/rnaseq_expression_classification.wdl), [BAM Workflow](./workflows/rnaseq_expression_classification/rnaseq_expression_classification_from_bams.wdl), [Counts Workflow](./workflows/rnaseq_expression_classification/rnaseq_expression_classification_from_counts.wdl) | ![In Production](https://img.shields.io/static/v1?label=Status&message=Production&color=green&style=flat-square)    |
| ESTIMATE                           | v1.0.0 (*beta*) | Runs the [ESTIMATE software package](https://bioinformatics.mdanderson.org/estimate/) on a feature counts file.                                       | None                                                                                                  | [Workflow](./workflows/rnaseq/ESTIMATE.wdl)                                                                                                                                                                                                                                                                                              | ![In Development](https://img.shields.io/static/v1?label=Status&message=Development&color=orange&style=flat-square) |
| Calculate Gene Lengths             | N/A             | Produces a gene length file from a GTF.                                                                                                               | None                                                                                                  | [Workflow](./workflows/rnaseq/calc-gene-lengths.wdl)                                                                                                                                                                                                                                                                                     | ![In Production](https://img.shields.io/static/v1?label=Status&message=Production&color=green&style=flat-square)    |
| Build BWA References               | N/A             | Builds reference files used by the [BWA aligner](https://github.com/lh3/bwa).                                                                         | None                                                                                                  | [Workflow](./workflows/general/bwa-db-build.wdl)                                                                                                                                                                                                                                                                                         | ![In Production](https://img.shields.io/static/v1?label=Status&message=Production&color=green&style=flat-square)    |
| BAM to FastQs                      | v1.0.0          | Split a BAM file into read groups, then read 1 FastQs and  read 2 FastQs.                                                                             | None                                                                                                  | [Workflow](./workflows/general/bam-to-fastqs.wdl)                                                                                                                                                                                                                                                                                        | ![In Production](https://img.shields.io/static/v1?label=Status&message=Production&color=green&style=flat-square)    |

## Author

👤 **St. Jude Cloud Team**

* Website: https://stjude.cloud
* Github: [@stjudecloud](https://github.com/stjudecloud)
* Twitter: [@StJudeResearch](https://twitter.com/StJudeResearch)

## Tests

Given that this repo is still new, there are no tests. When we add tests, we will update the README.

## 🤝 Contributing

Contributions, issues and feature requests are welcome!<br />Feel free to check [issues page](https://github.com/stjudecloud/workflows/issues). You can also take a look at the [contributing guide](https://github.com/stjudecloud/workflows/blob/master/CONTRIBUTING.md).

## Versioning

When versioned, workflows and Docker images will be versioned according to the [SemVer](http://semver.org/) guidelines. For now, we do not guarantee that all workflows will have an associated version.

## 📝 License

Copyright © 2019 [St. Jude Cloud Team](https://github.com/stjudecloud).<br />
This project is [MIT](https://github.com/stjudecloud/workflows/blob/master/LICENSE.md) licensed.

[wdl]: http://openwdl.org/
[cromwell]: https://github.com/broadinstitute/cromwell
[bioconda]: https://bioconda.github.io/
[oliver]: https://github.com/stjudecloud/oliver
