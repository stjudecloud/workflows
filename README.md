# St. Jude Cloud Workflows

This repository is intended to contain all source code related to bioinformatics
workflows in St. Jude Cloud. This includes but is not limited to Docker container
definition, WDL/CWL tools and workflows, and relevant documentation. The repo is 
currently under construction and should be considered in an experimental state. 

## RNA-seq workflow V2
The RNA-seq workflow used by St. Jude Cloud to generate BAM files is implemented 
in WDL. It uses the bioinformatics-base Docker image found in this repository. 
Configuration files for the Cromwell runner have been included to facilitate running 
on LSF in addition to local execution. 
