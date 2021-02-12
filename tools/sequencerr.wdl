

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

        mkdir ~{prefix}.sequencErr_results
        mv ~{prefix}.sequencErr.err ~{prefix}.sequencErr.count ~{prefix}.sequencErr_results
        tar -czf ~{prefix}.sequencErr_results.tar.gz ~{prefix}.sequencErr_results/
    }

    runtime {
        disk: disk_size + " GB"
        memory: "8 GB"
        docker: 'stjudecloud/sequencerr:branch-sequencErr-1.0.0'
        maxRetries: max_retries
    }

    output {
        File results = "~{prefix}.sequencErr_results.tar.gz"
    }
}