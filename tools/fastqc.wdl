## # FastQC
##
## This WDL tool wraps the [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) tool.
## FastQC generates quality control metrics for sequencing pipelines. 

version 1.0

task fastqc {
    input {
        File bam
        Int ncpu = 1
        Int memory_gb = 5
        Int max_retries = 1
    }

    String out_directory = basename(bam, ".bam") + ".fastqc_results"
    String out_tar_gz = out_directory + ".tar.gz"
    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil((bam_size * 2) + 10)

    command {
        set -euo pipefail
        
        mkdir ~{out_directory}
        fastqc -f bam \
            -o ~{out_directory} \
            -t ~{ncpu} \
            ~{bam}

        tar -czf ~{out_tar_gz} ~{out_directory}
    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        cpu: ncpu
        docker: 'ghcr.io/stjudecloud/fastqc:1.0.2'
        maxRetries: max_retries
    }

    output {
        File results = out_tar_gz
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "This WDL tool generates a FastQC quality control metrics report for the input BAM file."
    }

    parameter_meta {
        bam: "Input BAM format file to generate coverage for"
    }
}
