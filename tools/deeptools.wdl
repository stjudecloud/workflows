## # Deeptools
##
## This WDL tool wraps the [DeepTools](https://deeptools.readthedocs.io/en/develop/index.html) tool.
## DeepTools is a suite of Python tools for analysis of high throughput sequencing analysis.

version 1.0

task bamCoverage {
    input {
        File bam
        File bai
        String prefix = basename(bam, ".bam")
        Int max_retries = 1
        Int memory_gb = 5
        Int ncpu = 1
        Boolean detect_nproc = false
    }

    String parsed_detect_nproc = if detect_nproc then "true" else ""
    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil((bam_size * 4) + 10)
 
    command {
        set -euo pipefail
        
        n_cores=~{ncpu}
        if [ -n ~{parsed_detect_nproc} ]
        then
            n_cores="max"
        fi

        if [ ! -e ~{bam}.bai ]
        then 
            ln -s ~{bai} ~{bam}.bai
        fi
 
        bamCoverage --bam ~{bam} --outFileName ~{prefix}.bw --outFileFormat bigwig --numberOfProcessors "$n_cores"
    }

    runtime {
        cpu: ncpu
        disk: disk_size + " GB"
        memory: memory_gb + " GB"
        docker: 'stjudecloud/deeptools:1.0.1'
        maxRetries: max_retries
    }

    output {
        File bigwig = "~{prefix}.bw"
    }

    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
        description: "This WDL tool generates a BigWig coverage track using bamCoverage from DeepTools (https://deeptools.readthedocs.io/en/develop/index.html)."
    }

    parameter_meta {
        bam: "Input BAM format file to generate coverage for"
        bai: "BAM index file corresponding to the input BAM"
    }
} 
