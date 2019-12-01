# Workflows

This repository contains all bioinformatics workflows used on the St. Jude Cloud project. Officially, the repository is in beta — the project is adding workflows as they are developed and put into production.

## Getting Started

At the time of writing, all workflows are written in [WDL][wdl] and are tested using [Cromwell][cromwell]. Although we do not test outside of Cromwell, we expect that the workflows will work just as well using other runners.

The easiest way to get started is to install [bioconda][bioconda] and the run the following commands:

```bash
conda create -n workflows-dev -c bioconda cromwell==0.40 -y
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

* `conf` - Cromwell configuration files created for various environments that we use across our team. Feel free to use/fork/suggest improvements.
* `docker` - Dockerfiles used in our workflows. All docker images are published to [Docker Hub](https://hub.docker.com/u/stjudecloud).
  * At the time of writing, there is only one Docker image heavily in use, which is our [bioinformatics base image](./docker/bioinformatics-base/Dockerfile).
* `tools` - All tools we have wrapped as individual WDL tasks.
* `workflows` - Directory containing all end-to-end bioinformatics workflows.

## Workflows Available

The current workflows exist in this repo with the following statuses:

| Name             | Version         | Description                                            | Specification                                                                      | Workflow                                                    | Status                                                                                                              |
| ---------------- | --------------- | ------------------------------------------------------ | ---------------------------------------------------------------------------------- | ----------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------- |
| Build References | v2.0.0 *(beta)* | Build reference files used in harmonization pipelines. | None                                                                               | [Workflow](./workflows/reference/bootstrap-reference.wdl)   | ![In Development](https://img.shields.io/static/v1?label=Status&message=Development&color=orange&style=flat-square) |
| Standard RNA-Seq | v2.0.0 *(beta)* | Standard RNA-Seq harmonization pipeline.               | [Specification](https://stjudecloud.github.io/rfcs/0001-rnaseq-workflow-v2.0.html) | [Workflow](./workflows/standard-rnaseq/start-to-finish.wdl) | ![In Development](https://img.shields.io/static/v1?label=Status&message=Development&color=orange&style=flat-square) |

## Tests

Given that this repo is still new, there are no tests. When we add tests, we will update the README.

## Contributing

Please read [CONTRIBUTING.md](./CONTRIBUTING.md) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

When versioned, workflows will be versioned according to the [SemVer](http://semver.org/) guidelines. For now, we do not guarantee that all workflows will have an associated version.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

[wdl]: http://openwdl.org/
[cromwell]: https://github.com/broadinstitute/cromwell
[bioconda]: https://bioconda.github.io/