## # librarian

version 1.1

task librarian {
    meta {
        description: "Runs the `librarian` tool to derive the likely Illumina library preparation protocol used to generate a pair of FASTQ files."
        help: "**WARNING** this tool is not guaranteed to work on all data, and may produce nonsensical results. `librarian` was trained on a limited set of GEO read data (Gene Expression Oriented). This means the input data should be Paired-End, of mouse or human origin, read length should be >50bp, and derived from a library prep kit that is in the `librarian` database. This version of `librarian` has been trained on \"read one\" data of Paired-End sequencing data. It is not intended for use with Single-End data, even though it only accepts a single FASTQ."
        output: {
            report: "A tar archive containing the `librarian` report and raw data."
            raw_data: "The raw data that can be processed by MultiQC."
        }
    }

    parameter_meta {
        read_one_fastq: "Read one FASTQ of a Paired-End sample to analyze. May be uncompressed or gzipped."
        prefix: "Name of the output tar archive. The extension `.tar.gz` will be added."
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }
    input {
        File read_one_fastq
        String prefix = sub(
            basename(read_one_fastq),
            "([_\\.][rR][12])?(\\.subsampled)?\\.(fastq|fq)(\\.gz)?$",
            ""
        ) + ".librarian"
        Int modify_disk_size_gb = 0
    }

    Float read1_size = size(read_one_fastq, "GiB")
    Int disk_size_gb = (
        ceil(read1_size) + 10 + modify_disk_size_gb
    )

    command <<<
        set -euo pipefail

        mkdir ~{prefix}
        /app/librarian --local --raw -o ~{prefix} ~{read_one_fastq}

        tar -czf ~{prefix}.tar.gz ~{prefix}
    >>>

    output {
        File report = "~{prefix}.tar.gz"
        File raw_data = "~{prefix}/librarian_heatmap.txt"
    }

    runtime {
        memory: "4 GB"
        disks: "~{disk_size_gb} GB"
        docker: 'ghcr.io/desmondwillowbrook/librarian:1.2'
        maxRetries: 1
    }
}
