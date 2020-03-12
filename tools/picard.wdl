## Description:
##
## This WDL tool wraps the PicardTools library (https://broadinstitute.github.io/picard/).
## PicardTools is a set of Java tools for manipulating sequencing data.

version 1.0

task mark_duplicates {
    input {
        File bam
        String prefix = basename(bam, ".bam")
        Int memory_gb = 50
        Int max_retries = 1
        String wait_var = ""
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil((bam_size * 2) + 10)
    Int java_heap_size = ceil(memory_gb * 0.9)

    command {
        picard -Xmx~{java_heap_size}g MarkDuplicates I=~{bam} \
            O=~{prefix}.duplicates.bam \
            VALIDATION_STRINGENCY=SILENT \
            CREATE_INDEX=false \
            CREATE_MD5_FILE=false \
            COMPRESSION_LEVEL=5 \
            METRICS_FILE=~{prefix}.metrics.txt
    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        docker: 'stjudecloud/picard:1.0.0-alpha'
        maxRetries: max_retries
    }

    output {
        File out = "~{prefix}.duplicates.bam"
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "This WDL tool marks duplicate reads in the input BAM file using Picard."
    }

    parameter_meta {
        bam: "Input BAM format file to mark duplicates"
    }
}

task validate_bam {
    input {
        File input_bam
        Int memory_gb = 8
        Int max_retries = 1
        String wait_var = ""
    }

    Float bam_size = size(input_bam, "GiB")
    Int disk_size = ceil((bam_size * 2) + 10)
    Int java_heap_size = ceil(memory_gb * 0.9)
    
    command {
        picard -Xmx~{java_heap_size}g ValidateSamFile I=~{input_bam} \
            IGNORE=INVALID_PLATFORM_VALUE > stdout.txt
    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        docker: 'stjudecloud/picard:1.0.0-alpha'
        maxRetries: max_retries
    }

    output {
        String out = read_string("stdout.txt")
        File bam = input_bam
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "This WDL tool validates the input BAM file for correct formatting using Picard."
    }

    parameter_meta {
        bam: "Input BAM format file to validate"
    }
}

task bam_to_fastq {
    input {
        File bam
        String prefix = basename(bam, ".bam")
        Int memory_gb = 40
        Int max_retries = 1
        String wait_var = ""
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil((bam_size * 4) + 10)
    Int java_heap_size = ceil(memory_gb * 0.9)

    command {
        picard -Xmx~{java_heap_size}g SamToFastq INPUT=~{bam} \
            FASTQ=~{prefix}_R1.fastq \
            SECOND_END_FASTQ=~{prefix}_R2.fastq \
            RE_REVERSE=true \
            VALIDATION_STRINGENCY=SILENT
    }

    runtime{
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        docker: 'stjudecloud/picard:1.0.0-alpha'
        maxRetries: max_retries
    }

    output {
        File read1 = "~{prefix}_R1.fastq"
        File read2 = "~{prefix}_R2.fastq"
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "This WDL tool converts the input BAM file into paired FastQ format files."
    }

    parameter_meta {
        bam: "Input BAM format file to convert to FastQ"
    }
}

task sort {
    input {
        File bam
        String sort_order = "coordinate"
        String output_filename = basename(bam, ".bam") + ".sorted.bam"
        Int memory_gb = 25
        Int? disk_size_gb
        Int max_retries = 1
        String wait_var = ""
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = select_first([disk_size_gb, ceil((bam_size * 4) + 10)])
    Int java_heap_size = ceil(memory_gb * 0.9)

    command {
        picard -Xmx~{java_heap_size}g SortSam I=~{bam} \
           O=~{output_filename} \
           SO=~{sort_order} 
    }
    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        docker: 'stjudecloud/picard:1.0.0-alpha'
        maxRetries: max_retries
    }
    output {
        File sorted_bam = output_filename
    }
    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
        description: "This WDL tool sorts the input BAM file."
    }
    parameter_meta {
        bam: "Input BAM format file to sort"
    }
}
