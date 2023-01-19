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
        docker: 'ghcr.io/stjudecloud/picard:branch-replace_qualimap-1.1.0'
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
        Boolean succeed_on_errors = false
        Boolean succeed_on_warnings = true
        Array[String] ignore_list = ["MISSING_PLATFORM_VALUE", "INVALID_PLATFORM_VALUE", "INVALID_MAPPING_QUALITY"]
        Boolean summary_mode = false
        Boolean index_validation_stringency_less_exhaustive = false
        Int max_errors = 2147483647
        String output_filename = basename(bam, ".bam") + ".ValidateSamFile.txt"
        Int memory_gb = 8
        Int max_retries = 1
    }

    String succeed_on_errors_string = if (succeed_on_errors) then "true" else ""
    String succeed_on_warnings_string = if (succeed_on_warnings) then "true" else ""
    String mode_arg = if (summary_mode) then "MODE=SUMMARY" else ""
    String stringency_arg = if (index_validation_stringency_less_exhaustive)
        then "INDEX_VALIDATION_STRINGENCY=LESS_EXHAUSTIVE"
        else ""
    String ignore_prefix = if (length(ignore_list) != 0) then "IGNORE=" else ""

    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil((bam_size * 2) + 10)
    Int java_heap_size = ceil(memory_gb * 0.9)
    
    command {       
        picard -Xmx~{java_heap_size}g ValidateSamFile \
            I=~{bam} \
            ~{mode_arg} \
            ~{stringency_arg} \
            ~{ignore_prefix}~{sep=' IGNORE=' ignore_list} \
            MAX_OUTPUT=~{max_errors} \
            > ~{output_filename}

        rc=$?
        if [ $rc -le -1 ] || [ $rc -ge 4 ]; then
            exit $rc
        fi
        set -eo pipefail

        if [ "~{succeed_on_warnings_string}" == "true" ]; then
            GREP_PATTERN="ERROR"
        else
            GREP_PATTERN="(ERROR|WARNING)"
        fi

        if [ "~{succeed_on_errors_string}" != "true" ] && [ "$(grep -Ec "$GREP_PATTERN" ~{output_filename})" -gt 0 ]; then
            echo "Errors detected by Picard ValidateSamFile" > /dev/stderr
            grep -E "$GREP_PATTERN" ~{output_filename} > /dev/stderr
            exit 1
        fi
    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        docker: 'ghcr.io/stjudecloud/picard:branch-replace_qualimap-1.1.0'
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
        succeed_on_errors: "Succeed the task even if errors *and/or* warnings are detected"
        succeed_on_warnings: "Succeed the task if warnings are detected and there are no errors. Overridden by `succeed_on_errors`"
        ignore_list: "List of Picard errors and warnings to ignore. Possible values can be found on the [GATK website](https://gatk.broadinstitute.org/hc/en-us/articles/360035891231-Errors-in-SAM-or-BAM-files-can-be-diagnosed-with-ValidateSamFile)"
        summary_mode: "Boolean to enable SUMMARY mode of `picard ValidateSamFile`"
        index_validation_stringency_less_exhaustive: "Boolean to set `INDEX_VALIDATION_STRINGENCY=LESS_EXHAUSTIVE` for `picard ValidateSamFile`"
        max_errors: "Set the value of MAX_OUTPUT for `picard ValidateSamFile`. The Picard default is 100, a lower number can enable fast fail behavior"
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
        docker: 'ghcr.io/stjudecloud/picard:branch-replace_qualimap-1.1.0'
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
        docker: 'ghcr.io/stjudecloud/picard:branch-replace_qualimap-1.1.0'
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

task merge_sam_files {
    input {
        Array[File] bam
        String output_name = "merged.bam"
        String sort_order = "coordinate"
        Boolean threading = true
        Int memory_gb = 40
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil((bam_size * 2) + 10)
    Int java_heap_size = ceil(memory_gb * 0.9)
    Array[String] input_arg = prefix("INPUT=", bam)

    command {
        set -euo pipefail

        picard -Xmx~{java_heap_size}g MergeSamFiles \
            ~{sep=' ' input_arg} \
            OUTPUT=~{output_name} \
            SORT_ORDER=~{sort_order} \
            USE_THREADING=~{threading} \
            VALIDATION_STRINGENCY=SILENT
    }

    runtime{
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        docker: 'ghcr.io/stjudecloud/picard:branch-replace_qualimap-1.1.0'
        maxRetries: max_retries
    }

    output {
        File merged_bam = output_name
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "This WDL tool merges the input BAM files into a single BAM file."
    }

    parameter_meta {
        bam: "Input BAMs to merge"
    }
}

task clean_sam {
    input {
        File bam
        String output_filename = basename(bam, ".bam") + ".cleaned.bam"
        Int memory_gb = 25
        Int? disk_size_gb
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = select_first([disk_size_gb, ceil((bam_size * 2) + 10)])
    Int java_heap_size = ceil(memory_gb * 0.9)

    command {
        picard -Xmx~{java_heap_size}g CleanSam \
            I=~{bam} \
            O=~{output_filename}
    }
    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        docker: 'ghcr.io/stjudecloud/picard:branch-replace_qualimap-1.1.0'
        maxRetries: max_retries
    }
    output {
        File cleaned_bam = output_filename
    }
    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
        description: "This WDL tool cleans the input BAM file. Cleans soft-clipping beyond end-of-reference, sets MAPQ=0 for unmapped reads"
    }
    parameter_meta {
        bam: "Input BAM format file to clean"
    }
}