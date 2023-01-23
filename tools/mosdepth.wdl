## # Mosdepth
##
## This WDL tool wraps the [mosdepth tool](https://github.com/brentp/mosdepth).

version 1.0

task coverage {
    input {
        File bam
        File bai
        Int memory_gb = 8 
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil(bam_size + 5)

    command {
        set -euo pipefail

        mv ~{bai} ~{bam}.bai || true
        mosdepth -n -x "$(basename ~{bam} '.bam')" ~{bam}
    }

    runtime {
        disk: disk_size + " GB"
        memory: memory_gb + " GB"
        docker: 'ghcr.io/stjudecloud/mosdepth:branch-replace_qualimap-1.0.0'
        maxRetries: max_retries
    }

    output {
        File global_dist = basename(bam, '.bam') + ".mosdepth.global.dist.txt"
        File summary = basename(bam, '.bam') + ".mosdepth.summary.txt"
    }

    meta {
        author: "Andrew Frantz"
        email: "andrew.frantz@stjude.org"
        description: "This WDL tool wraps the Mosdepth tool for calculating coverage"
    }
}
