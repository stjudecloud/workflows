## [Homepage](https://lomereiter.github.io/sambamba/)
#
# SPDX-License-Identifier: MIT
# Copyright St. Jude Children's Research Hospital
version 1.1

task index {
    meta {
        description: "Creates a `.bai` BAM index for the input BAM"
        outputs: {
            bam_index: "A `.bai` BAM index associated with the input BAM. Filename will be `basename(bam) + '.bai'`."
        }
    }

    parameter_meta {
        bam: "Input BAM format file to index"
        use_all_cores: {
            description: "Use all cores? Recommended for cloud environments. Not recommended for cluster environments."
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
        Boolean use_all_cores = false
        Int ncpu = 2
        Int modify_disk_size_gb = 0
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size * 1.2) + 10 + modify_disk_size_gb

    String outfile_name = basename(bam) + ".bai"

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi

        sambamba index --nthreads "$n_cores" ~{bam} ~{outfile_name}
    >>>

    output {
        File bam_index = outfile_name
    }

    runtime {
        cpu: ncpu
        memory: "4 GB"
        disk: "~{disk_size_gb} GB"
        container: 'quay.io/biocontainers/sambamba:1.0--h98b6b92_0'
        maxRetries: 1
    }
}

task merge {
    meta {
        description: "Merges multiple sorted BAMs into a single BAM"
        outputs: {
            merged_bam: "The BAM resulting from merging all the input BAMs"
        }
    }

    parameter_meta {
        bams: "An array of BAMs to merge into one combined BAM"
        prefix: "Prefix for the BAM file. The extension `.bam` will be added."
        use_all_cores: {
            description: "Use all cores? Recommended for cloud environments. Not recommended for cluster environments."
            common: true
        }
        ncpu: {
            description: "Number of cores to allocate for task"
            common: true
        }
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        Array[File] bams
        String prefix
        Boolean use_all_cores = false
        Int ncpu = 2
        Int modify_disk_size_gb = 0
    }

    Float bams_size = size(bams, "GiB")
    Int disk_size_gb = ceil(bams_size * 2) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi

        sambamba merge \
            --nthreads "$n_cores" \
            ~{prefix}.bam \
            ~{sep(" ", bams)}
    >>>

    output {
        File merged_bam = prefix + ".bam"
    }

    runtime {
        cpu: ncpu
        memory: "4 GB"
        disk: "~{disk_size_gb} GB"
        container: 'quay.io/biocontainers/sambamba:1.0--h98b6b92_0'
        maxRetries: 1
    }
}

task sort {
    meta {
        description: "Sorts the input BAM file"
        outputs: {
            sorted_bam: "The input BAM after it has been sorted according to `sort_order`"
        }
    }

    parameter_meta {
        bam: "Input BAM format file to sort"
        sort_order: {
            description: "Order by which to sort the input BAM"
            choices: [
                'queryname',
                'coordinate'
            ]
            common: true
        }
        prefix: "Prefix for the sorted BAM file. The extension `.bam` will be added."
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        String sort_order = "coordinate"
        String prefix = basename(bam, ".bam") + ".sorted"
        Int memory_gb = 25
        Int modify_disk_size_gb = 0
        Int ncpu = 2
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size * 4) + 10 + modify_disk_size_gb

    String outfile_name = prefix + ".bam"

    command <<<
        set -euo pipefail

        sambamba sort \
            --nthreads ~{ncpu} \
            -o ~{outfile_name} \
            ~{if (sort_order == 'queryname') then '-n' else ''} \
            ~{bam}
    >>>

    output {
        File sorted_bam = outfile_name
    }

    runtime {
        cpu: ncpu
        memory: "~{memory_gb} GB"
        disk: "~{disk_size_gb} GB"
        container: 'quay.io/biocontainers/sambamba:1.0--h98b6b92_0'
        maxRetries: 1
    }
}

task markdup {
    meta {
        description: "Marks duplicate reads in the input BAM file using Picard"
        outputs: {
            duplicate_marked_bam: "The input BAM with computationally determined duplicates marked."
            mark_duplicates_metrics: "Duplicate marking metrics output from sambamba"
        }
    }

    parameter_meta {
        bam: "Input BAM format file in which to mark duplicates"
        prefix: "Prefix for the markdup result files. The extensions `.bam` will be added."
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        String prefix = basename(bam, ".bam") + ".MarkDuplicates"
        Boolean remove_duplicates = false
        Int memory_gb = 50
        Int modify_disk_size_gb = 0
        Int ncpu = 2
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil((bam_size * 2) + 10) + modify_disk_size_gb

    command <<<
        set -euo pipefail

        sambamba markdup \
            --nthreads ~{ncpu} \
            ~{if remove_duplicates then '--remove-duplicates' else ''} \
            ~{bam} \
            ~{prefix}.bam \
            > sambamba_markdup_log.txt
    >>>

    output {
        File duplicate_marked_bam = "~{prefix}.bam"
        File duplicate_marked_bam_index = "~{prefix}.bam.bai"
        File markdup_log = "sambamba_markdup_log.txt"
    }

    runtime {
        cpu: ncpu
        memory: "~{memory_gb} GB"
        disk: "~{disk_size_gb} GB"
        container: 'quay.io/biocontainers/sambamba:1.0--h98b6b92_0'
        maxRetries: 1
    }
}

task flagstat {
    meta {
        description: "Produces a samtools-like flagstat report containing statistics about the alignments based on the bit flags set in the BAM"
        outputs: {
            flagstat_report: "`sambamba flagstat` STDOUT redirected to a file"
        }
    }

    parameter_meta {
        bam: "Input BAM format file to generate flagstat for"
        outfile_name: "Name for the flagstat report file"
        use_all_cores: {
            description: "Use all cores? Recommended for cloud environments. Not recommended for cluster environments."
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
        String outfile_name = basename(bam, ".bam") + ".flagstat.txt"
        Boolean use_all_cores = false
        Int ncpu = 2
        Int modify_disk_size_gb = 0
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi

        sambamba flagstat --nthreads "$n_cores" ~{bam} > ~{outfile_name}
    >>>

    output {
       File flagstat_report = outfile_name
    }

    runtime {
        cpu: ncpu
        memory: "5 GB"
        disk: "~{disk_size_gb} GB"
        container: 'quay.io/biocontainers/sambamba:1.0--h98b6b92_0'
        maxRetries: 1
    }
}