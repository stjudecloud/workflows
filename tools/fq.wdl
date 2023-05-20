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
        docker: 'ghcr.io/stjudecloud/fqlib:1.2.0'
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

task subsample {
    input {
        File read1
        File? read2
        String prefix = basename(read1, "_R1.fastq.gz")
        Float probability = 1.0
        Int record_count = -1
        Int max_retries = 1
        Int memory_gb = 4
        Int modify_disk_size_gb = 0
    }

    parameter_meta {
        read1: "Input FastQ with read one"
        read2: "Input FastQ with read two"
        probability: "The probability a record is kept, as a percentage (0.0, 1.0). Cannot be used with `record-count`."
        record_count: "The exact number of records to keep. Cannot be used with `probability`."
        max_retries: "Number of times to retry in case of failure"
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    Float read1_size = size(read1, "GiB")
    Float read2_size = if defined(read2) then size(read2, "GiB") else 0

    Int disk_size_gb = ceil((read1_size + read2_size) * 2) + modify_disk_size_gb

    command <<<
        fq subsample \
            ~{if (probability < 1.0 && probability > 0)
                then "-p " + probability
                else ""
            } \
            ~{if (record_count > 0) then "-n " + record_count else ""} \
            --r1-dst ~{prefix + "_R1.subsampled.fastq.gz"} \
            ~{if defined(read2)
                then "--r2-dst " + prefix + "_R2.subsampled.fastq.gz"
                else ""
            } \
            ~{read1} \
            ~{if defined(read2) then read2 else ""}
    >>>

    runtime {
        disk: disk_size_gb + " GB"
        memory: memory_gb + " GB"
        docker: 'quay.io/biocontainers/fq:0.10.0--h9ee0642_0'
        maxRetries: max_retries
    }

    output {
        File subsampled_read1 = prefix + "_R1.subsampled.fastq.gz"
        File? subsampled_read2 = prefix + "_R2.subsampled.fastq.gz"
    }

    meta {
        author: "Andrew Frantz"
        email: "andrew.frantz@stjude.org"
        description: "This WDL task subsamples input FastQs"
    }
}
