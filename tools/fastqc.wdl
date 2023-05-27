## # FastQC
##
## This WDL file wraps the [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) tool.
## FastQC generates quality control metrics for sequencing pipelines. 

version 1.0

task fastqc {
    input {
        File bam
        Int ncpu = 1
        Int memory_gb = 5
        Boolean use_all_cores = false
        Int max_retries = 1
        Boolean use_all_cores = false
    }

    String out_directory = basename(bam, ".bam") + ".fastqc_results"
    String out_tar_gz = out_directory + ".tar.gz"

    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil((bam_size * 2) + 10)

    command {
        set -euo pipefail
        
        n_cores=~{ncpu}
        if [ "~{use_all_cores}" = "true" ]; then
            n_cores=$(nproc)
        fi
        
        mkdir ~{out_directory}
        fastqc -f bam \
            -o ~{out_directory} \
            -t "$n_cores" \
            ~{bam}

        tar -czf ~{out_tar_gz} ~{out_directory}
    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        cpu: ncpu
        docker: 'quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1'
        maxRetries: max_retries
    }

    output {
        File raw_data = "~{out_directory}/~{basename(bam, '.bam')}_fastqc.zip"
        File results = out_tar_gz
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "This WDL task generates a FastQC quality control metrics report for the input BAM file."
    }

    parameter_meta {
        bam: "Input BAM format file to generate coverage for"
    }
}
