# Building WDL workflows for DNAnexus

## Obtain dxWDL JAR
Retrieve the latest dxWDL JAR release from GitHub: https://github.com/dnanexus/dxWDL/releases

## Optional workflow parameters for dxWDL

* `-project` - Specify a project to compile the workflow. This is optional and otherwise uses the currently selected project.
* `-archive` - Archive older versions of the workflow and applets
* `-defaults` - Set default options for certain parameters
* `-verbose` - Detailed build information
* `-locked` - Creates a one stage worklfow that is cleaner in the interface
* `-extras` - JSON formatted file with options primarily for the DNAnexus platform settings

## Build Interactive t-SNE workflow for DNAnexus

Commands for building the t-SNE workflows are included below. Your version of dxWDL may differ from the version included below. Several optional parameters are included. `-defaults` specifies DNAnexus paths to reference data for the workflow. `-extras` specifies that tasks should be retried by default on failure. 

### Build workflow running htseq-count on BAM input

```
java -jar dxWDL-v1.46.2.jar compile workflows/interactive-tsne/interactive_tsne_from_bams.wdl -project project-FjFfvV89F80QvvxJ8131yzpB -archive -verbose -defaults workflows/interactive-tsne/inputs/defaults_bam.json -extras workflows/interactive-tsne/inputs/extras.json -locked
```

### Build workflow from HTSeq counts data

```
java -jar dxWDL-v1.46.2.jar compile workflows/interactive-tsne/interactive_tsne_from_counts.wdl -project project-FjFfvV89F80QvvxJ8131yzpB -archive -verbose -defaults workflows/interactive-tsne/inputs/defaults_counts.json -extras workflows/interactive-tsne/inputs/extras.json -locked
```

### Build workflow with RNA-Seq V2 remapping of BAM input

```
java -jar dxWDL-v1.46.2.jar compile workflows/interactive-tsne/interactive-tsne.wdl -project project-FjFfvV89F80QvvxJ8131yzpB -archive -verbose -defaults workflows/interactive-tsne/inputs/defaults.json  -extras workflows/interactive-tsne/inputs/extras.json -locked
```