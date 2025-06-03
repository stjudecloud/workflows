## [Homepage](https://lomereiter.github.io/sambamba/)

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
            description: "Use all cores? Recommended for cloud environments. Not recommended for cluster environments.",
            group: "Common",
        }
        ncpu: {
            description: "Number of cores to allocate for task",
            group: "Common",
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

        sambamba index --nthreads "$n_cores" "~{bam}" "~{outfile_name}"
    >>>

    output {
        File bam_index = outfile_name
    }

    runtime {
        cpu: ncpu
        memory: "4 GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/sambamba:1.0--h98b6b92_0"
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
            description: "Use all cores? Recommended for cloud environments. Not recommended for cluster environments.",
            group: "Common",
        }
        ncpu: {
            description: "Number of cores to allocate for task",
            group: "Common",
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
            "~{prefix}.bam" \
            ~{sep(" ", squote(bams))}
    >>>

    output {
        File merged_bam = prefix + ".bam"
    }

    runtime {
        cpu: ncpu
        memory: "4 GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/sambamba:1.0--h98b6b92_0"
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
        prefix: "Prefix for the sorted BAM file. The extension `.bam` will be added."
        queryname_sort: {
            description: "If true, sort the BAM by queryname. If false, sort by coordinate.",
            group: "Common",
        }
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
        ncpu: "Number of cores to allocate for task"
    }

    input {
        File bam
        String prefix = basename(bam, ".bam") + ".sorted"
        Boolean queryname_sort = false
        Int modify_disk_size_gb = 0
        Int ncpu = 2
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size * 3) + 10 + modify_disk_size_gb

    String outfile_name = prefix + ".bam"

    command <<<
        set -euo pipefail

        sambamba sort \
            --nthreads ~{ncpu} \
            -o "~{outfile_name}" \
            ~{if queryname_sort then "-n" else ""} \
            "~{bam}"
    >>>

    output {
        File sorted_bam = outfile_name
    }

    runtime {
        cpu: ncpu
        memory: "25 GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/sambamba:1.0--h98b6b92_0"
        maxRetries: 1
    }
}

task markdup {
    meta {
        description: "Marks duplicate reads in the input BAM file"
        outputs: {
            duplicate_marked_bam: "The input BAM with computationally determined duplicates marked.",
            duplicate_marked_bam_index: "Index file for the duplicate marked BAM",
            markdup_log: "Log file from the markdup process",
        }
    }

    parameter_meta {
        bam: "Input BAM format file in which to mark duplicates"
        prefix: "Prefix for the markdup result files. The extensions `markdup.bam` will be added."
        remove_duplicates: {
            description: "If true, remove duplicates instead of marking them.",
            group: "Common",
        }
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
        ncpu: "Number of cores to allocate for task"
    }

    input {
        File bam
        String prefix = basename(bam, ".bam")
        Boolean remove_duplicates = false
        Int modify_disk_size_gb = 0
        Int ncpu = 2
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil((bam_size * 2) + 10) + modify_disk_size_gb

    command <<<
        sambamba markdup \
            --nthreads ~{ncpu} \
            ~{if remove_duplicates then "--remove-duplicates" else ""} \
            "~{bam}" \
            "~{prefix}.markdup.bam" \
            > "~{prefix}.markdup_log.txt"
    >>>

    output {
        File duplicate_marked_bam = "~{prefix}.markdup.bam"
        File duplicate_marked_bam_index = "~{prefix}.markdup.bam.bai"
        File markdup_log = "~{prefix}.markdup_log.txt"
    }

    runtime {
        cpu: ncpu
        memory: "50 GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/sambamba:1.0--h98b6b92_0"
        maxRetries: 1
    }
}

task flagstat {
    meta {
        description: "Produces a report containing statistics about the alignments based on the bit flags set in the BAM"
        outputs: {
            flagstat_report: "`sambamba flagstat` STDOUT redirected to a file"
        }
    }

    parameter_meta {
        bam: "Input BAM format file to generate flagstat for"
        outfile_name: "Name for the flagstat report file"
        use_all_cores: {
            description: "Use all cores? Recommended for cloud environments. Not recommended for cluster environments.",
            group: "Common",
        }
        ncpu: {
            description: "Number of cores to allocate for task",
            group: "Common",
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

        sambamba flagstat --nthreads "$n_cores" "~{bam}" > "~{outfile_name}"
    >>>

    output {
       File flagstat_report = outfile_name
    }

    runtime {
        cpu: ncpu
        memory: "5 GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/sambamba:1.0--h98b6b92_0"
        maxRetries: 1
    }
}
