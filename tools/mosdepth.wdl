## [Homepage](https://github.com/brentp/mosdepth)
#
# SPDX-License-Identifier: MIT
# Copyright St. Jude Children's Research Hospital
version 1.1

task coverage {
    # TODO not all options exposed
    meta {
        description: "Runs the Mosdepth tool for calculating coverage"
        outputs: {
            summary: "A summary of mean depths per chromosome and within specified regions per chromosome"
            global_dist: "The `$prefix.mosdepth.global.dist.txt` file contains a cumulative distribution indicating the proportion of total bases that were covered for at least a given coverage value. It does this for each chromosome, and for the whole genome."
            region_dist: "The `$prefix.mosdepth.region.dist.txt` file contains a cumulative distribution indicating the proportion of total bases in the region(s) defined by `coverage_bed` that were covered for at least a given coverage value"
        }
    }

    parameter_meta {
        bam: "Input BAM format file to calculate coverage for"
        bam_index: "BAM index file corresponding to the input BAM"
        coverage_bed: "BED file to pass to the `-b` flag of `mosdepth`. This will restrict coverage analysis to regions defined by the BED file."
        prefix: "Prefix for the `mosdepth` report files. The extensions `.mosdepth.summary.txt`, `.mosdepth.global.dist.txt` and `.mosdepth.region.dist.txt` will be added."
        use_fast_mode: "Use Mosdepth's 'fast mode'? This enables the `-x` flag."
        min_mapping_quality: "Minimum mapping quality to pass to the `-Q` flag of `mosdepth`"
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
        max_retries: "Number of times to retry in case of failure"
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
        container: 'quay.io/biocontainers/mosdepth:0.3.3--h37c5b7d_2'
        maxRetries: max_retries
    }
}
