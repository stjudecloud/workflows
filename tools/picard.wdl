## # Picard
##
## This WDL tool wraps the [PicardTools library](https://broadinstitute.github.io/picard/).
## PicardTools is a set of Java tools for manipulating sequencing data.

version 1.0

task mark_duplicates {
    input {
        File bam
        String prefix = basename(bam, ".bam")
        Boolean create_bam = true
        Int memory_gb = 50
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = if create_bam then ceil((bam_size * 2) + 10) else ceil(bam_size + 10)
    Int java_heap_size = ceil(memory_gb * 0.9)

    command <<<
        set -euo pipefail

        picard -Xmx~{java_heap_size}g MarkDuplicates I=~{bam} \
            O=~{if create_bam then prefix + ".MarkDuplicates.bam" else "/dev/null"} \
            VALIDATION_STRINGENCY=SILENT \
            CREATE_INDEX=~{create_bam} \
            CREATE_MD5_FILE=~{create_bam} \
            COMPRESSION_LEVEL=5 \
            METRICS_FILE=~{prefix}.MarkDuplicates.metrics.txt
        
        mv ~{prefix}.MarkDuplicates.bai ~{prefix}.MarkDuplicates.bam.bai || true
    >>>

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        docker: 'ghcr.io/stjudecloud/picard:1.1.0'
        maxRetries: max_retries
    }

    output {
        File? duplicate_marked_bam = "~{prefix}.MarkDuplicates.bam"
        File? duplicate_marked_bam_index = "~{prefix}.MarkDuplicates.bam.bai"
        File? duplicate_marked_bam_md5 = "~{prefix}.MarkDuplicates.bam.md5"
        File mark_duplicates_metrics = "~{prefix}.MarkDuplicates.metrics.txt"
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
        String outfile_name = basename(bam, ".bam") + ".ValidateSamFile.txt"
        Int memory_gb = 12
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
            > ~{outfile_name}

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

        if [ "~{succeed_on_errors_string}" != "true" ] && [ "$(grep -Ec "$GREP_PATTERN" ~{outfile_name})" -gt 0 ]; then
            echo "Errors detected by Picard ValidateSamFile" > /dev/stderr
            grep -E "$GREP_PATTERN" ~{outfile_name} > /dev/stderr
            exit 1
        fi
    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        docker: 'ghcr.io/stjudecloud/picard:1.1.0'
        maxRetries: max_retries
    }

    output {
        File validate_report = outfile_name
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
        Boolean paired = true
        Int memory_gb = 56
        Int max_retries = 1
    }

    parameter_meta {
        bam: "Input BAM format file to convert to FastQ"
        paired: "Is the data paired-end (true) or single-end (false)?"
        max_retries: "Number of times to retry failed steps"
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil((bam_size * 4) + 10)
    Int java_heap_size = ceil(memory_gb * 0.9)

    command {
        set -euo pipefail

        picard -Xmx~{java_heap_size}g SamToFastq INPUT=~{bam} \
            FASTQ=~{prefix}_R1.fastq \
            ~{if paired then "SECOND_END_FASTQ="+prefix+"_R2.fastq" else ""} \
            RE_REVERSE=true \
            VALIDATION_STRINGENCY=SILENT
        
        gzip ~{prefix}_R1.fastq \
            ~{if paired then prefix+"_R2.fastq" else ""}
    }

    runtime{
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        docker: 'ghcr.io/stjudecloud/picard:1.1.0'
        maxRetries: max_retries
    }

    output {
        File read1 = "~{prefix}_R1.fastq.gz"
        File? read2 = "~{prefix}_R2.fastq.gz"
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "This WDL tool converts the input BAM file into paired FastQ format files."
    }
}

task sort {
    input {
        File bam
        String sort_order = "coordinate"
        String outfile_name = basename(bam, ".bam") + ".sorted.bam"
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
            O=~{outfile_name} \
            SO=~{sort_order} \
            CREATE_INDEX=false \
            CREATE_MD5_FILE=false \
            COMPRESSION_LEVEL=5 \
            VALIDATION_STRINGENCY=SILENT
    }
    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        docker: 'ghcr.io/stjudecloud/picard:1.1.0'
        maxRetries: max_retries
    }
    output {
        File sorted_bam = outfile_name
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
        String outfile_name = "merged.bam"
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
            OUTPUT=~{outfile_name} \
            SORT_ORDER=~{sort_order} \
            USE_THREADING=~{threading} \
            VALIDATION_STRINGENCY=SILENT
    }

    runtime{
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        docker: 'ghcr.io/stjudecloud/picard:1.1.0'
        maxRetries: max_retries
    }

    output {
        File merged_bam = outfile_name
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
        String outfile_name = basename(bam, ".bam") + ".cleaned.bam"
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
            O=~{outfile_name}
    }
    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        docker: 'ghcr.io/stjudecloud/picard:1.1.0'
        maxRetries: max_retries
    }
    output {
        File cleaned_bam = outfile_name
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

task collect_wgs_metrics {
    input {
        File bam
        File reference_fasta
        Int memory_gb = 12
        Int? disk_size_gb
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = select_first([disk_size_gb, ceil(bam_size + 5)])
    Int java_heap_size = ceil(memory_gb * 0.9)

    command {
        picard -Xmx~{java_heap_size}g CollectWgsMetrics \
            -I ~{bam} \
            -O "$(basename ~{bam} '.bam').CollectWgsMetrics.txt" \
            -R ~{reference_fasta} \
            --INCLUDE_BQ_HISTOGRAM true
    }
    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        docker: 'ghcr.io/stjudecloud/picard:1.1.0'
        maxRetries: max_retries
    }
    output {
        File wgs_metrics = basename(bam, ".bam") + ".CollectWgsMetrics.txt"
    }
    meta {
        author: "Andrew Frantz"
        email: "andrew.frantz@stjude.org"
        description: "This WDL tool runs the `picard CollectWgsMetrics` command."
    }
    parameter_meta {
        bam: "Input BAM format file to calculate WGS metrics for"
    }
}

task collect_wgs_metrics_with_nonzero_coverage {
    input {
        File bam
        File reference_fasta
        Int memory_gb = 12
        Int? disk_size_gb
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = select_first([disk_size_gb, ceil(bam_size + 5)])
    Int java_heap_size = ceil(memory_gb * 0.9)

    command {
        picard -Xmx~{java_heap_size}g CollectWgsMetricsWithNonZeroCoverage \
            -I ~{bam} \
            -O "$(basename ~{bam} '.bam').CollectWgsMetricsWithNonZeroCoverage.txt" \
            -CHART "$(basename ~{bam} '.bam').CollectWgsMetricsWithNonZeroCoverage.pdf" \
            -R ~{reference_fasta} \
            --INCLUDE_BQ_HISTOGRAM true
    }
    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        docker: 'ghcr.io/stjudecloud/picard:1.1.0'
        maxRetries: max_retries
    }
    output {
        File wgs_metrics = basename(bam, ".bam") + ".CollectWgsMetricsWithNonZeroCoverage.txt"
        File wgs_metrics_pdf = basename(bam, ".bam") + ".CollectWgsMetricsWithNonZeroCoverage.pdf"
    }
    meta {
        author: "Andrew Frantz"
        email: "andrew.frantz@stjude.org"
        description: "This WDL tool runs the `picard CollectWgsMetricsWithNonZeroCoverage` command."
    }
    parameter_meta {
        bam: "Input BAM format file to calculate WGS metrics for"
    }
}

task collect_alignment_summary_metrics {
    input {
        File bam
        Int memory_gb = 8
        Int? disk_size_gb
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = select_first([disk_size_gb, ceil(bam_size + 5)])
    Int java_heap_size = ceil(memory_gb * 0.9)

    command {
        picard -Xmx~{java_heap_size}g CollectAlignmentSummaryMetrics \
            -I ~{bam} \
            -O "$(basename ~{bam} '.bam').CollectAlignmentSummaryMetrics.txt" \
            -H "$(basename ~{bam} '.bam').CollectAlignmentSummaryMetrics.pdf"
    }
    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        docker: 'ghcr.io/stjudecloud/picard:1.1.0'
        maxRetries: max_retries
    }
    output {
        File alignment_metrics = basename(bam, ".bam") + ".CollectAlignmentSummaryMetrics.txt"
        File alignment_metrics_pdf = basename(bam, ".bam") + ".CollectAlignmentSummaryMetrics.pdf"
    }
    meta {
        author: "Andrew Frantz"
        email: "andrew.frantz@stjude.org"
        description: "This WDL tool runs the `picard CollectAlignmentSummaryMetrics` command."
    }
    parameter_meta {
        bam: "Input BAM format file to calculate alignment metrics for"
    }
}

task collect_gc_bias_metrics {
    input {
        File bam
        File reference_fasta
        Int memory_gb = 8
        Int? disk_size_gb
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = select_first([disk_size_gb, ceil(bam_size + 5)])
    Int java_heap_size = ceil(memory_gb * 0.9)

    command {
        picard -Xmx~{java_heap_size}g CollectGcBiasMetrics \
            -I ~{bam} \
            -R ~{reference_fasta} \
            -O "$(basename ~{bam} '.bam').CollectGcBiasMetrics.txt" \
            -S "$(basename ~{bam} '.bam').CollectGcBiasMetrics.summary.txt" \
            -CHART "$(basename ~{bam} '.bam').CollectGcBiasMetrics.pdf"
    }
    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        docker: 'ghcr.io/stjudecloud/picard:1.1.0'
        maxRetries: max_retries
    }
    output {
        File gc_bias_metrics = basename(bam, ".bam") + ".CollectGcBiasMetrics.txt"
        File gc_bias_metrics_pdf = basename(bam, ".bam") + ".CollectGcBiasMetrics.pdf"
    }
    meta {
        author: "Andrew Frantz"
        email: "andrew.frantz@stjude.org"
        description: "This WDL tool runs the `picard CollectGcBiasMetrics` command."
    }
    parameter_meta {
        bam: "Input BAM format file to calculate GC bias metrics for"
    }
}

task collect_insert_size_metrics {
    input {
        File bam
        Int memory_gb = 8
        Int? disk_size_gb
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = select_first([disk_size_gb, ceil(bam_size + 5)])
    Int java_heap_size = ceil(memory_gb * 0.9)

    command {
        picard -Xmx~{java_heap_size}g CollectInsertSizeMetrics \
            -I ~{bam} \
            -O "$(basename ~{bam} '.bam').CollectInsertSizeMetrics.txt" \
            -H "$(basename ~{bam} '.bam').CollectInsertSizeMetrics.pdf"
    }
    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        docker: 'ghcr.io/stjudecloud/picard:1.1.0'
        maxRetries: max_retries
    }
    output {
        File insert_size_metrics = basename(bam, ".bam") + ".CollectInsertSizeMetrics.txt"
        File insert_size_metrics_pdf = basename(bam, ".bam") + ".CollectInsertSizeMetrics.pdf"
    }
    meta {
        author: "Andrew Frantz"
        email: "andrew.frantz@stjude.org"
        description: "This WDL tool runs the `picard CollectInsertSizeMetrics` command."
    }
    parameter_meta {
        bam: "Input BAM format file to calculate insert size metrics for"
    }
}

task quality_score_distribution {
    input {
        File bam
        Int memory_gb = 8
        Int? disk_size_gb
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = select_first([disk_size_gb, ceil(bam_size + 5)])
    Int java_heap_size = ceil(memory_gb * 0.9)

    command {
        picard -Xmx~{java_heap_size}g QualityScoreDistribution \
            -I ~{bam} \
            -O "$(basename ~{bam} '.bam').QualityScoreDistribution.txt" \
            -CHART "$(basename ~{bam} '.bam').QualityScoreDistribution.pdf"
    }
    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        docker: 'ghcr.io/stjudecloud/picard:1.1.0'
        maxRetries: max_retries
    }
    output {
        File quality_score_distribution_txt = basename(bam, ".bam") + ".QualityScoreDistribution.txt"
        File quality_score_distribution_pdf = basename(bam, ".bam") + ".QualityScoreDistribution.pdf"
    }
    meta {
        author: "Andrew Frantz"
        email: "andrew.frantz@stjude.org"
        description: "This WDL tool runs the `picard QualityScoreDistribution` command."
    }
    parameter_meta {
        bam: "Input BAM format file to calculate quality score distribution for"
    }
}