## # Mosdepth
##
## This WDL file wraps the [mosdepth tool](https://github.com/brentp/mosdepth).

version 1.1

task coverage {
    meta {
        description: "This WDL task wraps the Mosdepth tool for calculating coverage"
    }

    input {
        File bam
        File bam_index
        File? coverage_bed
        String prefix = basename(bam, '.bam')
        Boolean use_fast_mode = true
        Int min_mapping_quality = 20
        Int memory_gb = 8 
        Int modify_disk_size_gb = 0
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        # localize BAM and BAI to CWD
        # some backends prevent writing to the inputs directories
        # to accomodate this, create symlinks in CWD
        CWD_BAM=~{basename(bam)}
        ln -s ~{bam} "$CWD_BAM"
        ln -s ~{bam_index} "$CWD_BAM".bai

        mosdepth \
            -n \
            ~{if defined(coverage_bed) then "-b" else ""} ~{coverage_bed} \
            -Q ~{min_mapping_quality} \
            ~{if (use_fast_mode) then "-x" else ""} \
            ~{prefix} \
            "$CWD_BAM"

        rm "$CWD_BAM" "$CWD_BAM".bai
    >>>

    output {
        File summary = prefix + ".mosdepth.summary.txt"
        File global_dist = prefix + ".mosdepth.global.dist.txt"
        File? region_dist = prefix + ".mosdepth.region.dist.txt"
    }

    runtime {
        memory: "~{memory_gb} GB"
        disk: "~{disk_size_gb} GB"
        docker: 'quay.io/biocontainers/mosdepth:0.3.3--h37c5b7d_2'
        maxRetries: max_retries
    }
}
