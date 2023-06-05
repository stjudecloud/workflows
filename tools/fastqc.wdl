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
        String prefix = basename(bam, ".bam") + ".fastqc_results"
        Boolean use_all_cores = false
        Int ncpu = 1
        Int memory_gb = 5
        Int modify_disk_size_gb = 0
        Int max_retries = 1
    }

    String out_tar_gz = prefix + ".tar.gz"

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size * 2) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail
        
        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi
        
        mkdir ~{prefix}
        fastqc -f bam \
            -o ~{prefix} \
            -t "$n_cores" \
            ~{bam}

        tar -czf ~{out_tar_gz} ~{prefix}
    >>>

    output {
        File raw_data = "~{prefix}/~{basename(bam, '.bam')}_fastqc.zip"  # TODO verify this works if prefix differs
        File results = out_tar_gz
    }

    runtime {
        cpu: ncpu
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        docker: 'quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1'
        maxRetries: max_retries
    }
}
