

version 1.0

task sequencerr {
    input {
        File bam
        File bai
        String prefix = basename(bam, ".bam")
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil(bam_size + 8)

    command {
        set -euo pipefail
        
        mv ~{bai} ~{bam}.bai || true
        sequencerr -pe=~{prefix}.sequencErr.err ~{bam} ~{prefix}.sequencErr.count
    }

    runtime {
        disk: disk_size + " GB"
        memory: "8 GB"
        docker: 'stjudecloud/branch-sequencErr:1.0.0'
        maxRetries: max_retries
    }

    output {
        File error = "~{prefix}.sequencErr.err"
        File count = "~{prefix}.sequencErr.count"
    }
}