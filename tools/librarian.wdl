## # librarian

task librarian {
    input {
        File read_one_fastq
        File read_two_fastq
        String prefix = sub(
            basename(read_one_fastq),
            "([_\.][rR][12])?(\.subsampled)?\.(fastq|fq)(\.gz)?$",
            ""
        )
        Int memory_gb = 4
        Int modify_disk_size_gb = 0
        Int max_retries = 1
    }

    Float read1_size = size(read_one_fastq, "GiB")
    Float read2_size = size(read_two_fastq, "GiB")
    Int disk_size_gb = (
        ceil(read1_size + read2_size) + 10 + modify_disk_size_gb
    )

    command <<<
        /app/librarian --local -o ~{prefix} ~{read_one_fastq} ~{read_two_fastq}
    >>>

    # output {
    #     File report = 
    # }

    runtime {
        memory: "~{memory_gb} GB"
        disk: "~{disk_size_gb} GB"
        docker: 'ghcr.io/desmondwillowbrook/librarian:latest'
        maxRetries: max_retries
    }
}
