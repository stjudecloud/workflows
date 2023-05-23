## # FQ
##
## This WDL file wraps the [fq tool](https://github.com/stjude/fqlib).
## The fq library provides methods for manipulating Illumina generated 
## FastQ files.

version 1.0

task fqlint {
    input {
        File read1
        File? read2
        Int max_retries = 1
        Int modify_memory_gb = 0
    }

    Float read1_size = size(read1, "GiB")
    Float read2_size = if defined(read2) then size(read2, "GiB") else 0

    Int memory_gb_calculation = ceil(((read1_size + read2_size) * 0.08)) + modify_memory_gb
    Int memory_gb = if memory_gb_calculation > 4
        then memory_gb_calculation
        else 4

    Int disk_size = ceil((read1_size + read2_size) * 2)

    String args = if defined(read2) then "" else "--disable-validator P001" 

    command {
        fq lint ~{args} ~{read1} ~{read2}
    }

    runtime {
        disk: disk_size + " GB"
        memory: memory_gb + " GB"
        docker: 'quay.io/biocontainers/fq:0.9.1--h9ee0642_0'
        maxRetries: max_retries
    }

    output {
        File validated_read1 = read1
        File? validated_read2 = read2
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "This WDL task performs quality control on the input FastQ pairs to ensure proper formatting."
    }

    parameter_meta {
        read1: "Input FastQ with read one"
        read2: "Input FastQ with read two"
    }
}
