## # FQ
##
## This WDL tool wraps the [fq tool](https://github.com/stjude/fqlib).
## The fq library provides methods for manipulating Illumina generated 
## FastQ files.

version 1.0

task fqlint {
    input {
        File read1
        File? read2
        Int max_retries = 1
        Int memory_gb = 8 
    }

    Float read1_size = size(read1, "GiB")
    Float read2_size = if defined(read2) then size(read2, "GiB") else 0
    Int disk_size = ceil(((read1_size + read2_size) * 2) + 10)
    String args = if defined(read2) then "" else "--disable-validator P001" 

    command {
        fq lint ~{args} ~{read1} ~{read2}
    }

    runtime {
        disk: disk_size + " GB"
        memory: memory_gb + " GB"
        docker: 'ghcr.io/stjudecloud/fqlib:1.0.1'
        maxRetries: max_retries
    }

    output {
        File validated_read1 = read1
        File? validated_read2 = read2
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "This WDL tool performs quality control on the input FastQ pairs to ensure proper formatting."
    }

    parameter_meta {
        read1: "Input FastQ with read one"
        read2: "Input FastQ with read two"
    }
}
