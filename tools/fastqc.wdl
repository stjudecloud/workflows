## # FastQC
##
## This WDL file wraps the [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) tool.
## FastQC generates quality control metrics for sequencing pipelines. 

version 1.0

task fastqc {
    meta {
        description: "This WDL task generates a FastQC quality control metrics report for the input BAM file."
    }

    parameter_meta {
        bam: "Input BAM format file to run FastQC on"
    }

    input {
        File bam
        String results_directory = basename(bam, ".bam") + ".fastqc_results"
        Boolean use_all_cores = false
        Int memory_gb = 5
        Int modify_disk_size_gb = 0
        Int ncpu = 1
        Int max_retries = 1
    }

    String out_tar_gz = results_directory + ".tar.gz"

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil((bam_size * 2) + 10) + modify_disk_size_gb

    command <<<
        set -euo pipefail
        
        n_cores=~{ncpu}
        if [ "~{use_all_cores}" = "true" ]; then
            n_cores=$(nproc)
        fi
        
        mkdir ~{results_directory}
        fastqc -f bam \
            -o ~{results_directory} \
            -t "$n_cores" \
            ~{bam}

        tar -czf ~{out_tar_gz} ~{results_directory}
    >>>

    output {
        File raw_data = "~{results_directory}/~{basename(bam, '.bam')}_fastqc.zip"
        File results = out_tar_gz
    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        cpu: ncpu
        docker: 'quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1'
        maxRetries: max_retries
    }
}
