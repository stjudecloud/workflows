## Description:
##
## This WDL tool wraps the DeepTools tool (https://deeptools.readthedocs.io/en/develop/index.html).
## DeepTools is a suite of Python tools for analysis of high throughput sequencing analysis.

version 1.0

task bamCoverage {
    input {
        File bam
        File bai
        String prefix = basename(bam, ".bam")
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil((bam_size * 4) + 10)
 
    command {
        if [ ! -e ${bam}.bai ] 
        then 
            ln -s ${bai} ${bam}.bai
        fi
 
        bamCoverage --bam ${bam} --outFileName ${prefix}.bw --outFileFormat bigwig --numberOfProcessors "max"
    }

    runtime {
        disk: disk_size + " GB"
        docker: 'stjudecloud/deeptools:1.0.0-alpha'
        maxRetries: max_retries
    }

    output {
        File bigwig = "${prefix}.bw"
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
