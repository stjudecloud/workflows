## # Mosdepth
##
## This WDL file wraps the [mosdepth tool](https://github.com/brentp/mosdepth).

version 1.0

task coverage {
    input {
        File bam
        File bai
        Int min_mapping_quality = 20
        Boolean use_fast_mode = true
        Int memory_gb = 8 
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil(bam_size + 5)

    command {
        set -euo pipefail

        mv ~{bai} ~{bam}.bai || true
        mosdepth \
            -n \
            -Q ~{min_mapping_quality} \
            ~{if (use_fast_mode) then "-x" else ""} \
            "$(basename ~{bam} '.bam')" \
            ~{bam}
    }

    runtime {
        disk: disk_size + " GB"
        memory: memory_gb + " GB"
        docker: 'quay.io/biocontainers/mosdepth:0.3.3--h37c5b7d_2'
        maxRetries: max_retries
    }

    output {
        File global_dist = basename(bam, '.bam') + ".mosdepth.global.dist.txt"
        File summary = basename(bam, '.bam') + ".mosdepth.summary.txt"
    }

    meta {
        author: "Andrew Frantz"
        email: "andrew.frantz@stjude.org"
        description: "This WDL task wraps the Mosdepth tool for calculating coverage"
    }
}
