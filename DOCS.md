<p align="center">
  <a href="https://github.com/stjudecloud/workflows"><img src="./assets/workflows-banner-flowchart.jpg" width="800" title="St. Jude Cloud Workflows"></a>
  <br />
  <a href="https://github.com/stjudecloud/workflows/blob/master/LICENSE.md" target="_blank">
    <img alt="License: MIT" src="https://img.shields.io/badge/License-MIT-yellow.svg" />
  </a>
</p>

> This repository contains all bioinformatics workflows used on the St. Jude Cloud project. Officially, the repository is in beta — the project is adding workflows as they are developed and put into production.

Resources requirements have been optimized to minimize failures in our computing environment, but they may not reflect the best settings for your use case. Please ensure that you tailor these parameters to fit your needs.

### 🏠 [Homepage](https://stjude.cloud)

## Workflows Structure

The repository contains workflows from multiple categories:

* `Harmonization` - These workflows are standard entrypoints used by the St. Jude Cloud team to harmonize data for release on stjude.cloud.
* `Reference` - These workflows are used to generate the necessary reference file inputs for our other workflows.
* `Utility` - These workflows are generally called as subworkflows of other workflows, but may have utility on their own.
* `Variant Calling` - These workflows are used to call a set of variants for a sample.
* `Other` - These workflows subworkflows of other workflows that should not be called in isolation or experimental workflows.
* `External` - These workflows are hosted in other repositories, but are used by one (or more) of the workflows in this repository.

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

👤 **St. Jude Cloud Team**

* Website: https://stjude.cloud
* Github: [@stjudecloud](https://github.com/stjudecloud)
* Twitter: [@StJudeResearch](https://twitter.com/StJudeResearch)

## 🤝 Contributing

Contributions, issues and feature requests are welcome!<br />Feel free to check [issues page](https://github.com/stjudecloud/workflows/issues). You can also take a look at the [contributing guide](https://github.com/stjudecloud/workflows/blob/master/CONTRIBUTING.md).

## 📝 License

Copyright © 2020-Present [St. Jude Cloud Team](https://github.com/stjudecloud).<br />
This project is [MIT](https://github.com/stjudecloud/workflows/blob/master/LICENSE.md) licensed.
