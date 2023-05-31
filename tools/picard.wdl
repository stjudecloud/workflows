## # Picard
##
## This WDL file wraps the [PicardTools library](https://broadinstitute.github.io/picard/).
## PicardTools is a set of Java tools for manipulating sequencing data.

version 1.0

task mark_duplicates {
    meta {
        description: "This WDL task marks duplicate reads in the input BAM file using Picard."
    }

    parameter_meta {
        bam: "Input BAM format file to mark duplicates"
    }

    input {
        File bam
        String prefix = basename(bam, ".bam") + ".MarkDuplicates"
        Boolean create_bam = true
        Int memory_gb = 50
        Int modify_disk_size_gb = 0
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = (
        (
            if create_bam
            then ceil((bam_size * 2) + 10)
            else ceil(bam_size + 10)
        ) + modify_disk_size_gb
    )
    
    Int java_heap_size = ceil(memory_gb * 0.9)

    command <<<
        set -euo pipefail

        picard -Xmx~{java_heap_size}g MarkDuplicates I=~{bam} \
            O=~{if create_bam then prefix + ".bam" else "/dev/null"} \
            VALIDATION_STRINGENCY=SILENT \
            CREATE_INDEX=~{create_bam} \
            CREATE_MD5_FILE=~{create_bam} \
            COMPRESSION_LEVEL=5 \
            METRICS_FILE=~{prefix}.metrics.txt
        
        if ~{create_bam}; then
            mv ~{prefix}.bai ~{prefix}.bam.bai
        fi
    >>>

    output {
        File? duplicate_marked_bam = "~{prefix}.bam"
        File? duplicate_marked_bam_index = "~{prefix}.bam.bai"
        File? duplicate_marked_bam_md5 = "~{prefix}.bam.md5"
        File mark_duplicates_metrics = "~{prefix}.metrics.txt"
    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        docker: 'quay.io/biocontainers/picard:2.27.5--hdfd78af_0'
        maxRetries: max_retries
    }
}

task validate_bam {
    # TODO should this be refactored to behave as "default" Picard behaves?
    #   Default Picard has some weird/not ideal behaviors
    #   e.g. `max_errors = 100`
    meta {
        description: "This WDL task validates the input BAM file for correct formatting using Picard."
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

    input {
        File bam
        Array[String] ignore_list = ["MISSING_PLATFORM_VALUE", "INVALID_PLATFORM_VALUE", "INVALID_MAPPING_QUALITY"]
        String outfile_name = basename(bam, ".bam") + ".ValidateSamFile.txt"
        Boolean succeed_on_errors = false
        Boolean succeed_on_warnings = true
        Boolean summary_mode = false
        Boolean index_validation_stringency_less_exhaustive = false
        Int max_errors = 2147483647  # max 32-bit INT
        Int memory_gb = 16
        Int modify_disk_size_gb = 0
        Int max_retries = 1
    }

    String mode_arg = if (summary_mode) then "MODE=SUMMARY" else ""
    String stringency_arg = if (index_validation_stringency_less_exhaustive)
        then "INDEX_VALIDATION_STRINGENCY=LESS_EXHAUSTIVE"
        else ""
    String ignore_prefix = if (length(ignore_list) != 0) then "IGNORE=" else ""

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size * 2) + 10 + modify_disk_size_gb
    Int java_heap_size = ceil(memory_gb * 0.9)
    
    command <<<
        picard -Xmx~{java_heap_size}g ValidateSamFile \
            I=~{bam} \
            ~{mode_arg} \
            ~{stringency_arg} \
            ~{ignore_prefix}~{sep=' IGNORE=' ignore_list} \
            MAX_OUTPUT=~{max_errors} \
            > ~{outfile_name}

        rc=$?
        if [ $rc -le -1 ] || [ $rc -ge 4 ]; then  # TODO explain this
            exit $rc
        fi

        set -euo pipefail

        if ~{succeed_on_warnings}; then
            GREP_PATTERN="ERROR"
        else
            GREP_PATTERN="(ERROR|WARNING)"
        fi

        if ~{succeed_on_errors} \
            && [ "$(grep -Ec "$GREP_PATTERN" ~{outfile_name})" -gt 0 ]
        then
            echo "Errors detected by Picard ValidateSamFile" > /dev/stderr
            grep -E "$GREP_PATTERN" ~{outfile_name} > /dev/stderr
            exit 1
        fi
    >>>

    output {
        File validate_report = outfile_name
        File validated_bam = bam
    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        docker: 'quay.io/biocontainers/picard:2.27.5--hdfd78af_0'
        maxRetries: max_retries
    }
}

task bam_to_fastq {
    meta {
        description: "*[Deprecated]* This WDL task converts the input BAM file into FastQ format files. This task has been deprecated in favor of `samtools.collate_to_fastq` which is more performant and doesn't error on 'illegal mate states'."
    }

    parameter_meta {
        bam: "Input BAM format file to convert to FastQ"
        paired: "Is the data paired-end (true) or single-end (false)?"
        max_retries: "Number of times to retry failed steps"
    }

    input {
        File bam
        String prefix = basename(bam, ".bam")
        Boolean paired = true
        Int memory_gb = 56
        Int modify_disk_size_gb = 0
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size * 4) + 10 + modify_disk_size_gb
    Int java_heap_size = ceil(memory_gb * 0.9)

    command <<<
        set -euo pipefail

        picard -Xmx~{java_heap_size}g SamToFastq INPUT=~{bam} \
            FASTQ=~{prefix}_R1.fastq \
            ~{if paired then "SECOND_END_FASTQ="+prefix+"_R2.fastq" else ""} \
            RE_REVERSE=true \
            VALIDATION_STRINGENCY=SILENT
        
        gzip ~{prefix}_R1.fastq \
            ~{if paired then prefix+"_R2.fastq" else ""}
    >>>

    output {
        File read_one_fastq_gz = "~{prefix}_R1.fastq.gz"
        File? read_two_fastq_gz = "~{prefix}_R2.fastq.gz"
    }

    runtime{
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        docker: 'quay.io/biocontainers/picard:2.27.5--hdfd78af_0'
        maxRetries: max_retries
    }
}

task sort {
    meta {
        description: "This WDL task sorts the input BAM file."
    }

    parameter_meta {
        bam: "Input BAM format file to sort"
    }

    input {
        File bam
        String sort_order = "coordinate"
        String prefix = basename(bam, ".bam") + ".sorted"
        Int memory_gb = 25
        Int modify_disk_size_gb = 0
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size * 4) + 10 + modify_disk_size_gb
    Int java_heap_size = ceil(memory_gb * 0.9)

    String outfile_name = prefix + ".bam"

    command <<<
        picard -Xmx~{java_heap_size}g SortSam \
            I=~{bam} \
            O=~{outfile_name} \
            SO=~{sort_order} \
            CREATE_INDEX=false \
            CREATE_MD5_FILE=false \
            COMPRESSION_LEVEL=5 \
            VALIDATION_STRINGENCY=SILENT
    >>>

    output {
        File sorted_bam = outfile_name
    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        docker: 'quay.io/biocontainers/picard:2.27.5--hdfd78af_0'
        maxRetries: max_retries
    }
}

task merge_sam_files {
    meta {
        description: "This WDL task merges the input BAM files into a single BAM file."
    }

    parameter_meta {
        bam: "Input BAMs to merge"
        threading: "Option to create a background thread to encode, compress and write to disk the output file. The threaded version uses about 20% more CPU and decreases runtime by ~20% when writing out a compressed BAM file."
    }

    input {
        Array[File] bams
        String prefix
        String sort_order = "coordinate"
        Boolean threading = true
        Int memory_gb = 40
        Int modify_disk_size_gb = 0
        Int max_retries = 1
    }

    Float bams_size = size(bams, "GiB")
    Int disk_size_gb = ceil(bams_size * 2) + 10 + modify_disk_size_gb
    Int java_heap_size = ceil(memory_gb * 0.9)

    Array[String] input_arg = prefix("INPUT=", bams)

    String outfile_name = prefix + ".bam"

    command <<<
        picard -Xmx~{java_heap_size}g MergeSamFiles \
            ~{sep=' ' input_arg} \
            OUTPUT=~{outfile_name} \
            SORT_ORDER=~{sort_order} \
            USE_THREADING=~{threading} \
            VALIDATION_STRINGENCY=SILENT
    >>>

    runtime{
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        docker: 'quay.io/biocontainers/picard:2.27.5--hdfd78af_0'
        maxRetries: max_retries
    }

    output {
        File merged_bam = outfile_name
    }
}

task clean_sam {
    meta {
        description: "This WDL task cleans the input BAM file. Cleans soft-clipping beyond end-of-reference, sets MAPQ=0 for unmapped reads"
    }

    parameter_meta {
        bam: "Input BAM format file to clean"
    }

    input {
        File bam
        String prefix = basename(bam, ".bam") + ".cleaned"
        Int memory_gb = 25
        Int modify_disk_size_gb = 0
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size * 2) + 10 + modify_disk_size_gb
    Int java_heap_size = ceil(memory_gb * 0.9)

    String outfile_name = prefix + ".bam"

    command <<<
        picard -Xmx~{java_heap_size}g CleanSam \
            I=~{bam} \
            O=~{outfile_name}
    >>>

    output {
        File cleaned_bam = outfile_name
    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        docker: 'quay.io/biocontainers/picard:2.27.5--hdfd78af_0'
        maxRetries: max_retries
    }
}

task collect_wgs_metrics {
    meta {
        description: "This WDL task runs the `picard CollectWgsMetrics` command."
    }

    parameter_meta {
        bam: "Input BAM format file to calculate WGS metrics for"
    }

    input {
        File bam
        File reference_fasta
        String outfile_name = basename(bam, ".bam") + ".CollectWgsMetrics.txt"
        Int memory_gb = 12
        Int modify_disk_size_gb = 0
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size) + 10 + modify_disk_size_gb
    Int java_heap_size = ceil(memory_gb * 0.9)

    command <<<
        picard -Xmx~{java_heap_size}g CollectWgsMetrics \
            -I ~{bam} \
            -R ~{reference_fasta} \
            -O ~{outfile_name} \
            --INCLUDE_BQ_HISTOGRAM true 
    >>>

    output {
        File wgs_metrics = outfile_name
    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        docker: 'quay.io/biocontainers/picard:2.27.5--hdfd78af_0'
        maxRetries: max_retries
    }
}

task collect_wgs_metrics_with_nonzero_coverage {
    meta {
        description: "*[Deprecated]* This WDL task runs the `picard CollectWgsMetricsWithNonZeroCoverage` command. This command is labelled as 'EXPERIMENTAL' in the Picard documentation."  # TODO does more need to be said?
    }

    parameter_meta {
        bam: "Input BAM format file to calculate WGS metrics for"
    }

    input {
        File bam
        File reference_fasta
        String prefix = basename(bam, ".bam")
        Int memory_gb = 12
        Int modify_disk_size_gb = 0
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size) + 10 + modify_disk_size_gb
    Int java_heap_size = ceil(memory_gb * 0.9)

    command <<<
        picard -Xmx~{java_heap_size}g CollectWgsMetricsWithNonZeroCoverage \
            -I ~{bam} \
            -O ~{prefix}.CollectWgsMetricsWithNonZeroCoverage.txt \
            -CHART ~{prefix}.CollectWgsMetricsWithNonZeroCoverage.pdf \
            -R ~{reference_fasta} \
            --INCLUDE_BQ_HISTOGRAM true
    >>>

    output {
        File wgs_metrics = prefix + ".CollectWgsMetricsWithNonZeroCoverage.txt"
        File wgs_metrics_pdf = prefix + ".CollectWgsMetricsWithNonZeroCoverage.pdf"
    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        docker: 'quay.io/biocontainers/picard:2.27.5--hdfd78af_0'
        maxRetries: max_retries
    }
}

task collect_alignment_summary_metrics {
    meta {
        author: "Andrew Frantz"
        email: "andrew.frantz@stjude.org"
        description: "This WDL task runs the `picard CollectAlignmentSummaryMetrics` command."
    }

    parameter_meta {
        bam: "Input BAM format file to calculate alignment metrics for"
    }

    input {
        File bam
        String prefix = basename(bam, ".bam") + ".CollectAlignmentSummaryMetrics"
        Int memory_gb = 8
        Int modify_disk_size_gb = 0
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size) + 10 + modify_disk_size_gb
    Int java_heap_size = ceil(memory_gb * 0.9)

    command <<<
        picard -Xmx~{java_heap_size}g CollectAlignmentSummaryMetrics \
            -I ~{bam} \
            -O ~{prefix}.txt \
            -H ~{prefix}.pdf
    >>>

    output {
        File alignment_metrics = prefix + ".txt"
        File alignment_metrics_pdf = prefix + ".pdf"
    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        docker: 'quay.io/biocontainers/picard:2.27.5--hdfd78af_0'
        maxRetries: max_retries
    }
}

task collect_gc_bias_metrics {
    meta {
        description: "This WDL task runs the `picard CollectGcBiasMetrics` command."
    }

    parameter_meta {
        bam: "Input BAM format file to calculate GC bias metrics for"
    }

    input {
        File bam
        File reference_fasta
        String prefix = basename(bam, ".bam") + ".CollectGcBiasMetrics"
        Int memory_gb = 8
        Int modify_disk_size_gb = 0
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size) + 10 + modify_disk_size_gb
    Int java_heap_size = ceil(memory_gb * 0.9)

    command <<<
        picard -Xmx~{java_heap_size}g CollectGcBiasMetrics \
            -I ~{bam} \
            -R ~{reference_fasta} \
            -O ~{prefix}.txt \
            -S ~{prefix}.summary.txt \
            -CHART ~{prefix}.pdf
    >>>

    output {
        File gc_bias_metrics = prefix + ".txt"
        File gc_bias_metrics_summary = prefix + ".summary.txt"
        File gc_bias_metrics_pdf = prefix + ".pdf"
    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        docker: 'quay.io/biocontainers/picard:2.27.5--hdfd78af_0'
        maxRetries: max_retries
    }
}

task collect_insert_size_metrics {
    meta {
        description: "This WDL task runs the `picard CollectInsertSizeMetrics` command."
    }

    parameter_meta {
        bam: "Input BAM format file to calculate insert size metrics for"
    }

    input {
        File bam
        String prefix = basename(bam, ".bam") + ".CollectInsertSizeMetrics"
        Int memory_gb = 8
        Int modify_disk_size_gb = 0
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size) + 10 + modify_disk_size_gb
    Int java_heap_size = ceil(memory_gb * 0.9)

    command <<<
        picard -Xmx~{java_heap_size}g CollectInsertSizeMetrics \
            -I ~{bam} \
            -O ~{prefix}.txt \
            -H ~{prefix}.pdf
    >>>

    output {
        File insert_size_metrics = prefix + ".txt"
        File insert_size_metrics_pdf = prefix + ".pdf"
    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        docker: 'quay.io/biocontainers/picard:2.27.5--hdfd78af_0'
        maxRetries: max_retries
    }
}

task quality_score_distribution {
    meta {
        description: "This WDL task runs the `picard QualityScoreDistribution` command."
    }

    parameter_meta {
        bam: "Input BAM format file to calculate quality score distribution for"
    }

    input {
        File bam
        String prefix = basename(bam, ".bam") + ".QualityScoreDistribution"
        Int memory_gb = 8
        Int modify_disk_size_gb = 0
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size) + 10 + modify_disk_size_gb
    Int java_heap_size = ceil(memory_gb * 0.9)

    command <<<
        picard -Xmx~{java_heap_size}g QualityScoreDistribution \
            -I ~{bam} \
            -O ~{prefix}.txt \
            -CHART ~{prefix}.pdf
    >>>

    output {
        File quality_score_distribution_txt = prefix + ".txt"
        File quality_score_distribution_pdf = prefix + ".pdf"
    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        docker: 'quay.io/biocontainers/picard:2.27.5--hdfd78af_0'
        maxRetries: max_retries
    }
}
