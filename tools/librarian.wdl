## # librarian

version 1.1

task librarian {
    input {
        File read_one_fastq
        String prefix = sub(
            basename(read_one_fastq),
            "([_\\.][rR][12])?(\\.subsampled)?\\.(fastq|fq)(\\.gz)?$",
            ""
        )
        Int modify_disk_size_gb = 0
        Int max_retries = 1
    }

    Float read1_size = size(read_one_fastq, "GiB")
    Int disk_size_gb = (
        ceil(read1_size) + 10 + modify_disk_size_gb
    )

    command <<<
        set -euo pipefail

        mkdir tmp
        TMPDIR=$(pwd)/tmp RUST_LOG=trace /app/librarian --local --raw -o ~{prefix} ~{read_one_fastq}
    >>>

    # output {
    #     File report =
    # }

    runtime {
        memory: "4 GB"
        disk: "~{disk_size_gb} GB"
        docker: 'ghcr.io/desmondwillowbrook/librarian:1.2'
        maxRetries: max_retries
    }
}
