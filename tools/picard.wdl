## # Picard
##
## This WDL tool wraps the [PicardTools library](https://broadinstitute.github.io/picard/).
## PicardTools is a set of Java tools for manipulating sequencing data.

version 1.0

task mark_duplicates {
    input {
        File bam
        String prefix = basename(bam, ".bam")
        Int memory_gb = 50
        Int max_retries = 1
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
        docker: 'stjudecloud/picard:1.0.1'
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
        File bam
        Boolean pass_on_errors = false
        Boolean pass_on_warnings = true
        Array[String] ignore_list = ["MISSING_PLATFORM_VALUE", "INVALID_PLATFORM_VALUE", "INVALID_MAPPING_QUALITY"]
        Boolean summary_mode = false
        Boolean index_validation_stringency_less_exhaustive = false
        Int max_errors = 100
        String output_filename = basename(bam, ".bam") + ".ValidateSamFile.txt"
        Int memory_gb = 8
        Int max_retries = 1
    }

    String pass_on_errors_string = if (pass_on_errors) then "true" else ""
    String pass_on_warnings_string = if (pass_on_warnings) then "true" else ""
    String mode_arg = if (summary_mode) then "MODE=SUMMARY" else ""
    String stringency_arg = if (index_validation_stringency_less_exhaustive)
        then "INDEX_VALIDATION_STRINGENCY=LESS_EXHAUSTIVE"
        else ""

    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil((bam_size * 2) + 10)
    Int java_heap_size = ceil(memory_gb * 0.9)
    
    command {       
        set -eo pipefail

        picard -Xmx~{java_heap_size}g ValidateSamFile \
            I=~{bam} \
            ~{mode_arg} \
            ~{stringency_arg} \
            ~{sep=' IGNORE=' ignore_list} \
            MAX_OUTPUT=~{max_errors} \
            > ~{output_filename}

        if [ "~{pass_on_warnings_string}" == "true" ]; then
            GREP_PATTERN="ERROR"
        else
            GREP_PATTERN="(ERROR|WARNING)"
        fi

        if [ "~{pass_on_errors_string}" != "true" ] && [ "$(grep -Ec "$GREP_PATTERN" ~{output_filename})" -gt 0 ]; then
            echo "Errors detected by Picard ValidateSamFile" > /dev/stderr
            grep -E "$GREP_PATTERN" ~{output_filename} > /dev/stderr
            exit 1
        fi
    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        docker: 'stjudecloud/picard:1.0.1'
        maxRetries: max_retries
    }

    output {
        File out = output_filename
        File validated_bam = bam
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
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil((bam_size * 4) + 10)
    Int java_heap_size = ceil(memory_gb * 0.9)

    command {
        set -euo pipefail

        picard -Xmx~{java_heap_size}g SamToFastq INPUT=~{bam} \
            FASTQ=~{prefix}_R1.fastq \
            SECOND_END_FASTQ=~{prefix}_R2.fastq \
            RE_REVERSE=true \
            VALIDATION_STRINGENCY=SILENT
        
        gzip ~{prefix}_R1.fastq ~{prefix}_R2.fastq
    }

    runtime{
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        docker: 'stjudecloud/picard:1.0.1'
        maxRetries: max_retries
    }

    output {
        File read1 = "~{prefix}_R1.fastq.gz"
        File read2 = "~{prefix}_R2.fastq.gz"
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
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = select_first([disk_size_gb, ceil((bam_size * 4) + 10)])
    Int java_heap_size = ceil(memory_gb * 0.9)

    command {
        picard -Xmx~{java_heap_size}g SortSam \
            I=~{bam} \
            O=~{output_filename} \
            SO=~{sort_order} \
            CREATE_INDEX=false \
            CREATE_MD5_FILE=false \
            COMPRESSION_LEVEL=5 \
            VALIDATION_STRINGENCY=SILENT
    }
    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        docker: 'stjudecloud/picard:1.0.1'
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
