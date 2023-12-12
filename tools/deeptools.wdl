## [Homepage](https://deeptools.readthedocs.io/en/develop/index.html)
#
# SPDX-License-Identifier: MIT
# Copyright St. Jude Children's Research Hospital
version 1.1

task bam_coverage {
    # TODO expose other params/formats
    meta {
        description: "Generates a BigWig coverage track using bamCoverage from DeepTools"
        outputs: {
            bigwig: "BigWig format coverage file"
        }
    }

    parameter_meta {
        bam: "Input BAM format file to generate coverage for"
        bam_index: "BAM index file corresponding to the input BAM"
        prefix: "Prefix for the BigWig file. The extension `.bw` will be added."
        use_all_cores: {
            description: "Use all cores? Recommended for cloud environments."
            common: true
        }
        ncpu: {
            description: "Number of cores to allocate for task"
            common: true
        }
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        File bam_index
        String prefix = basename(bam, ".bam")
        Boolean use_all_cores = false
        Int ncpu = 4
        Int modify_disk_size_gb = 0
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size * 1.5) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores="max"
        fi

        # localize BAM and BAI to CWD
        CWD_BAM=~{basename(bam)}
        ln -s ~{bam} "$CWD_BAM"
        ln -s ~{bam_index} "$CWD_BAM".bai

        bamCoverage \
            --bam "$CWD_BAM" \
            --outFileName ~{prefix}.bw \
            --outFileFormat bigwig \
            --numberOfProcessors "$n_cores"

        rm "$CWD_BAM" "$CWD_BAM".bai
    >>>

    output {
        File bigwig = "~{prefix}.bw"
    }

    runtime {
        cpu: ncpu
        memory: "4 GB"
        disk: "~{disk_size_gb} GB"
        container: 'quay.io/biocontainers/deeptools:3.5.1--pyhdfd78af_1'
        maxRetries: 1
    }
}
