<p align="center">
  <a href="https://github.com/stjudecloud/workflows"><img src="./docs/workflows-banner-flowchart.jpg" width="800" title="St. Jude Cloud Workflows"></a>
  <br />
  <a href="https://actions-badge.atrox.dev/stjudecloud/workflows/goto?ref=main">
    <img alt="Build Status" src="https://img.shields.io/endpoint.svg?url=https%3A%2F%2Factions-badge.atrox.dev%2Fstjudecloud%2Fworkflows%2Fbadge%3Fref%3Dmain&style=flat" />
  </a>
  <a href="https://stjudecloud.github.io/workflows/" target="_blank">
    <img alt="Documentation" src="https://img.shields.io/badge/documentation-yes-brightgreen.svg" />
  </a>
  <a href="https://github.com/stjudecloud/workflows/blob/master/LICENSE.md" target="_blank">
    <img alt="License: MIT" src="https://img.shields.io/badge/License-MIT-yellow.svg" />
  </a>
</p>

> This repository contains all bioinformatics workflows used on the St. Jude Cloud project. Officially, the repository is in beta — the project is adding workflows as they are developed and put into production.

Resources requirements have been optimized to minimize failures in our computing environment, but they may not reflect the best settings for your use case. Please ensure that you tailor these parameters to fit your needs.

### 🏠 [Homepage](https://stjude.cloud)

Please excuse the state of our documentation. We are working on some big changes around here, and with those changes will come much improved documentation.

## Repository Structure

The repository is laid out as follows:

* `workflows/` - Directory containing all end-to-end bioinformatics workflows.
* `tools/` - All tools we have wrapped as individual WDL tasks.
* `docker/` - Dockerfiles used in our workflows. All docker images are published to the [GitHub Container Registy](https://github.com/orgs/stjudecloud/packages?repo_name=workflows) as a part of our CI and are versioned.
* `tests/tools/` - Home to all of our testing infrastructure. We use [pytest-workflow](https://pytest-workflow.readthedocs.io/en/stable/) for validating our code.
* `bin/` - **no longer in use** Scripts used by Cromwell configuration settings. Add this to `$PATH` prior to using configurations in `conf` with Cromwell.
* `conf/` - **no longer in use** Cromwell configuration files created for various environments that we use across our team. Feel free to use/fork/suggest improvements.

## Author

👤 **St. Jude Cloud Team**

* Website: https://stjude.cloud
* Github: [@stjudecloud](https://github.com/stjudecloud)
* Twitter: [@StJudeResearch](https://twitter.com/StJudeResearch)

## Tests

Every task in this repository is covered by at least one test (see all of our tests in `tests/tools/`). These are run using [pytest-workflow](https://pytest-workflow.readthedocs.io/en/stable/).

## 🤝 Contributing

Contributions, issues and feature requests are welcome!<br />Feel free to check [issues page](https://github.com/stjudecloud/workflows/issues). You can also take a look at the [contributing guide](https://github.com/stjudecloud/workflows/blob/master/CONTRIBUTING.md).

## 📝 License

Copyright © 2020-Present [St. Jude Cloud Team](https://github.com/stjudecloud).<br />
This project is [MIT](https://github.com/stjudecloud/workflows/blob/master/LICENSE.md) licensed.
