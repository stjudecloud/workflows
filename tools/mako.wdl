version 1.3

enum SortOrder[String] {
    queryname = "queryname",
    coordinate = "coordinate",
    queryname_natural = "queryname::natural",
    template_coordinate = "template-coordinate",
}

task sort {
    meta {
        description: "Sorts the input BAM file"
        outputs: {
            sorted_bam: "The input BAM after it has been sorted according to `sort_order`",
            sorted_bam_index: "The index file for the sorted BAM file, if `write_index` is true",
        }
    }

    parameter_meta {
        bam: "Input BAM format file to sort"
        sort_order: {
            description: "Order by which to sort the input BAM",
            choices: [
                "queryname",
                "coordinate",
                "queryname::natural",
                "template-coordinate",
            ],
            group: "Common",
        }
        prefix: "Prefix for the sorted BAM file and accessory files. The extension `.bam` will be added."
        write_index: "Write a BAM index file for the sorted BAM file."
        verify: "Only verify sort order. Does not sort the BAM file."
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
        ncpu: "Number of CPUs to allocate for task"
    }

    input {
        File bam
        String sort_order = "coordinate"
        String prefix = basename(bam, ".bam") + ".sorted"
        Boolean write_index = true
        Boolean verify = false
        Int memory_gb = 25
        Int modify_disk_size_gb = 0
        Int ncpu = 1
    }

    Float bam_size = size(bam, "GB")
    Int disk_size_gb = ceil(bam_size * 4) + 10 + modify_disk_size_gb

    String outfile_name = prefix + ".bam"

    Int task_mem_gb = max(memory_gb - 2, 2)  # Reserve 2GB for overhead

    command <<<
        set -euo pipefail

        mako sort \
            ~{if verify
                then "--verify"
                else "-o \"~{outfile_name}\""
            } \
            ~{if write_index
                then "--write-index"
                else ""
            } \
            --order "~{sort_order}" \
            --max-memory "~{task_mem_gb}GB" \
            --threads "~{ncpu}" \
            -i "~{bam}"
    >>>

    output {
        File sorted_bam = outfile_name
        File? sorted_bam_index = outfile_name + ".bai"
    }

    requirements {
        cpu: ncpu
        memory: "~{memory_gb} GB"
        disks: "~{disk_size_gb} GB"
        container: "ghcr.io/stjudecloud/mako:0.1.3-0"
        maxRetries: 1
    }
}
