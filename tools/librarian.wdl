## # librarian

version 1.1

task librarian {
    meta {
        description: "Runs the Librarian tool to derive the likely Illumina library preparation protocol used to generate a read one FASTQ file."
        help: "This version of Librarian has been trained on \"read one\" data of Paired-End sequencing data. It is not intended for use with Single-End data, even though it only accepts a single FASTQ."
        output: {
            File report = "A tar archive containing the Librarian report and raw data."
            File raw_data = "The raw data that can be passed to MultiQC."
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
        RUST_LOG=trace /app/librarian --local --raw -o ~{prefix} ~{read_one_fastq}

        tar -czf ~{prefix}.tar.gz ~{prefix}
    >>>

    output {
        File report = "~{prefix}.tar.gz"
        File raw_data = "~{prefix}/librarian_heatmap.txt"
    }

    runtime {
        memory: "4 GB"
        disk: "~{disk_size_gb} GB"
        docker: 'ghcr.io/desmondwillowbrook/librarian:1.2'
        maxRetries: 1
    }
}
