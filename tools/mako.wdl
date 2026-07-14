version 1.4

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
        verify: "Only verify sort order. Does not sort the BAM file."
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        String sort_order = "coordinate"
        String prefix = basename(bam, ".bam") + ".sorted"
        Boolean verify = false
        Int memory_gb = 25
        Int modify_disk_size_gb = 0
    }

    Float bam_size = size(bam, "GB")
    Int disk_size_gb = ceil(bam_size * 4) + 10 + modify_disk_size_gb

    String outfile_name = prefix + ".bam"

    command <<<
        set -euo pipefail

        mako sort \
            -i "~{bam}" \
            ~{if verify then "--verify" else "-o \"~{outfile_name}\""} \
            --order "~{sort_order}"
    >>>

    output {
        File sorted_bam = outfile_name
    }

    requirements {
        memory: "~{memory_gb} GB"
        disks: "~{disk_size_gb} GB"
        container: "ghcr.io/stjudecloud/mako:0.1.3-0"
        maxRetries: 1
    }
}