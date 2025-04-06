## [Homepage](https://broadinstitute.github.io/picard/)

version 1.1

task validate_bam {
    meta {
        description: "Validates the input BAM file for correct formatting using Picard"
        external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360057440611-ValidateSamFile-Picard-"
        outputs: {
            validate_report: "Validation report produced by `picard ValidateSamFile`. Validation warnings and errors are logged.",
       }
   }

    parameter_meta {
        bam: "Input BAM format file to validate"
        reference_fasta: "Reference genome in FASTA format. Presence of the reference FASTA allows for `NM` tag validation."
        ignore_list: {
            description: "List of Picard errors and warnings to ignore. Possible values can be found on the GATK website (see `external_help`).",
            external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360035891231-Errors-in-SAM-or-BAM-files-can-be-diagnosed-with-ValidateSamFile",
            group: "common",
       }
        outfile_name: "Name for the ValidateSamFile report file"
        validation_stringency: {
            description: "Validation stringency for parsing the input BAM.",
            choices: [
                "STRICT",
                "LENIENT",
                "SILENT"
            ],
            tool_default: "STRICT",
       }
        succeed_on_errors: {
            description: "Succeed the task even if errors *and/or* warnings are detected",
            group: "common",
       }
        succeed_on_warnings: {
            description: "Succeed the task if warnings are detected and there are no errors. Overridden by `succeed_on_errors`",
            group: "common",
       }
        summary_mode: {
            description: "Enable SUMMARY mode?",
            group: "common",
       }
        index_validation_stringency_less_exhaustive: "Set `INDEX_VALIDATION_STRINGENCY=LESS_EXHAUSTIVE`?"
        max_errors: "Set the value of MAX_OUTPUT for `picard ValidateSamFile`. The Picard default is 100, a lower number can enable fast fail behavior"
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
   }

    input {
        File bam
        File? reference_fasta
        Array[String] ignore_list = []
        String outfile_name = basename(bam, ".bam") + ".ValidateSamFile.txt"
        String validation_stringency = "LENIENT"
        Boolean succeed_on_errors = false
        Boolean succeed_on_warnings = true
        Boolean summary_mode = false
        Boolean index_validation_stringency_less_exhaustive = false
        Int max_errors = 2147483647  # max 32-bit INT
        Int memory_gb = 16
        Int modify_disk_size_gb = 0
   }

    String reference_arg = (
        if defined(reference_fasta)
        then "-R ~{reference_fasta}"
        else ""
    )
    String mode_arg = if (summary_mode) then "--MODE SUMMARY" else ""
    String stringency_arg = (
        if (index_validation_stringency_less_exhaustive)
        then "--INDEX_VALIDATION_STRINGENCY LESS_EXHAUSTIVE"
        else ""
    )
    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size * 2) + 10 + modify_disk_size_gb
    Int java_heap_size = ceil(memory_gb * 0.9)

    command <<<
        set -euo pipefail

        rc=0
        picard -Xmx~{java_heap_size}g ValidateSamFile \
            -I ~{bam} \
            ~{reference_arg} \
            ~{mode_arg} \
            ~{stringency_arg} \
            --VALIDATION_STRINGENCY ~{validation_stringency} \
            ~{sep(" ", prefix("--IGNORE ", ignore_list))} \
            --MAX_OUTPUT ~{max_errors} \
            > ~{outfile_name} \
            || rc=$?

        # rc = 0 = success
        # rc = 1 = validation warnings (no errors)
        # rc = 2 = validation errors and warnings
        # rc = 3 = validation errors (no warnings)
        if [ $rc -ne 0 ] && [ $rc -ne 1 ] && [ $rc -ne 2 ] && [ $rc -ne 3 ]; then
            exit $rc
        fi

        if ~{succeed_on_warnings}; then
            GREP_PATTERN="ERROR"
        else
            GREP_PATTERN="(ERROR|WARNING)"
        fi

        if ! ~{succeed_on_errors} \
            && [ "$(grep -Ec "$GREP_PATTERN" ~{outfile_name})" -gt 0 ]
        then
            >&2 echo "Problems detected by Picard ValidateSamFile"
            >&2 grep -E "$GREP_PATTERN" ~{outfile_name}
            exit $rc
        fi
    >>>

    output {
        File validate_report = outfile_name
   }

    runtime {
        memory: "~{memory_gb} GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/picard:3.1.1--hdfd78af_0"
        maxRetries: 1
   }
}

task collect_wgs_metrics {
    meta {
        description: "Runs `picard CollectWgsMetrics` to collect metrics about the fractions of reads that pass base- and mapping-quality filters as well as coverage (read-depth) levels"
        external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360037226132-CollectWgsMetrics-Picard-"
        outputs: {
            wgs_metrics: {
                description: "Output report of `picard CollectWgsMetrics`",
                external_help: "https://broadinstitute.github.io/picard/picard-metric-definitions.html#CollectWgsMetrics.WgsMetrics",
           }
       }
   }

    parameter_meta {
        bam: "Input BAM format file for which to calculate WGS metrics"
        reference_fasta: "Gzipped reference genome in FASTA format"
        outfile_name: "Name for the metrics result file"
        validation_stringency: {
            description: "Validation stringency for parsing the input BAM.",
            choices: [
                "STRICT",
                "LENIENT",
                "SILENT"
            ],
            tool_default: "STRICT",
       }
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
   }

    input {
        File bam
        File reference_fasta
        String outfile_name = basename(bam, ".bam") + ".CollectWgsMetrics.txt"
        String validation_stringency = "SILENT"
        Int memory_gb = 12
        Int modify_disk_size_gb = 0
   }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size) + 10 + modify_disk_size_gb
    Int java_heap_size = ceil(memory_gb * 0.9)

    command <<<
        picard -Xmx~{java_heap_size}g CollectWgsMetrics \
            -I ~{bam} \
            -R ~{reference_fasta} \
            -O ~{outfile_name} \
            --VALIDATION_STRINGENCY ~{validation_stringency} \
            --INCLUDE_BQ_HISTOGRAM true
    >>>

    output {
        File wgs_metrics = outfile_name
   }

    runtime {
        memory: "~{memory_gb} GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/picard:3.1.1--hdfd78af_0"
        maxRetries: 1
   }
}

task collect_alignment_summary_metrics {
    meta {
        description: "Runs `picard CollectAlignmentSummaryMetrics` to calculate metrics detailing the quality of the read alignments as well as the proportion of the reads that passed machine signal-to-noise threshold quality filters"
        external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360040507751-CollectAlignmentSummaryMetrics-Picard-"
        outputs: {
            alignment_metrics: {
                description: "The text file output of `CollectAlignmentSummaryMetrics`",
                external_help: "http://broadinstitute.github.io/picard/picard-metric-definitions.html#AlignmentSummaryMetrics",
            },
            alignment_metrics_pdf: "The PDF file output of `CollectAlignmentSummaryMetrics`",
       }
   }

    parameter_meta {
        bam: "Input BAM format file for which to calculate alignment metrics"
        prefix: "Prefix for the output report files. The extensions `.txt` and `.pdf` will be added."
        validation_stringency: {
            description: "Validation stringency for parsing the input BAM.",
            choices: [
                "STRICT",
                "LENIENT",
                "SILENT"
            ],
            tool_default: "STRICT",
       }
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
   }

    input {
        File bam
        String prefix = basename(bam, ".bam") + ".CollectAlignmentSummaryMetrics"
        String validation_stringency = "SILENT"
        Int memory_gb = 8
        Int modify_disk_size_gb = 0
   }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size) + 10 + modify_disk_size_gb
    Int java_heap_size = ceil(memory_gb * 0.9)

    command <<<
        picard -Xmx~{java_heap_size}g CollectAlignmentSummaryMetrics \
            -I ~{bam} \
            --VALIDATION_STRINGENCY ~{validation_stringency} \
            -O ~{prefix}.txt \
            -H ~{prefix}.pdf
    >>>

    output {
        File alignment_metrics = prefix + ".txt"
        File alignment_metrics_pdf = prefix + ".pdf"
   }

    runtime {
        memory: "~{memory_gb} GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/picard:3.1.1--hdfd78af_0"
        maxRetries: 1
   }
}

task collect_gc_bias_metrics {
    meta {
        description: "Runs `picard CollectGcBiasMetrics` to collect information about the relative proportions of guanine (G) and cytosine (C) nucleotides"
        external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360037593931-CollectGcBiasMetrics-Picard-"
        outputs: {
            gc_bias_metrics: {
                description: "The full text file output of `CollectGcBiasMetrics`",
                external_help: "http://broadinstitute.github.io/picard/picard-metric-definitions.html#GcBiasDetailMetrics",
            },
            gc_bias_metrics_summary: {
                description: "The summary text file output of `CollectGcBiasMetrics`",
                external_help: "http://broadinstitute.github.io/picard/picard-metric-definitions.html#GcBiasSummaryMetrics",
            },
            gc_bias_metrics_pdf: "The PDF file output of `CollectGcBiasMetrics`",
       }
   }

    parameter_meta {
        bam: "Input BAM format file for which to calculate GC bias metrics"
        reference_fasta: "Reference sequences in FASTA format"
        prefix: "Prefix for the output report files. The extensions `.txt`, `.summary.txt`, and `.pdf` will be added."
        validation_stringency: {
            description: "Validation stringency for parsing the input BAM.",
            choices: [
                "STRICT",
                "LENIENT",
                "SILENT"
            ],
            tool_default: "STRICT",
       }
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
   }

    input {
        File bam
        File reference_fasta
        String prefix = basename(bam, ".bam") + ".CollectGcBiasMetrics"
        String validation_stringency = "SILENT"
        Int memory_gb = 8
        Int modify_disk_size_gb = 0
   }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size) + 10 + modify_disk_size_gb
    Int java_heap_size = ceil(memory_gb * 0.9)

    command <<<
        picard -Xmx~{java_heap_size}g CollectGcBiasMetrics \
            -I ~{bam} \
            -R ~{reference_fasta} \
            --VALIDATION_STRINGENCY ~{validation_stringency} \
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
        memory: "~{memory_gb} GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/picard:3.1.1--hdfd78af_0"
        maxRetries: 1
   }
}

task collect_insert_size_metrics {
    meta {
        description: "Runs `picard CollectInsertSizeMetrics` to collect metrics for validating library construction including the insert size distribution and read orientation of Paired-End libraries"
        external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360037055772-CollectInsertSizeMetrics-Picard-"
        outputs: {
            insert_size_metrics: {
                description: "The text file output of `CollectInsertSizeMetrics`",
                external_help: "http://broadinstitute.github.io/picard/picard-metric-definitions.html#InsertSizeMetrics",
            },
            insert_size_metrics_pdf: "The PDF file output of `CollectInsertSizeMetrics`",
       }
   }

    parameter_meta {
        bam: "Input BAM format file for which to calculate insert size metrics"
        prefix: "Prefix for the output report files. The extensions `.txt` and `.pdf` will be added."
        validation_stringency: {
            description: "Validation stringency for parsing the input BAM.",
            choices: [
                "STRICT",
                "LENIENT",
                "SILENT"
            ],
            tool_default: "STRICT",
       }
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
   }

    input {
        File bam
        String prefix = basename(bam, ".bam") + ".CollectInsertSizeMetrics"
        String validation_stringency = "SILENT"
        Int memory_gb = 8
        Int modify_disk_size_gb = 0
   }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size) + 10 + modify_disk_size_gb
    Int java_heap_size = ceil(memory_gb * 0.9)

    command <<<
        picard -Xmx~{java_heap_size}g CollectInsertSizeMetrics \
            -I ~{bam} \
            --VALIDATION_STRINGENCY ~{validation_stringency} \
            -O ~{prefix}.txt \
            -H ~{prefix}.pdf
    >>>

    output {
        File insert_size_metrics = prefix + ".txt"
        File insert_size_metrics_pdf = prefix + ".pdf"
   }

    runtime {
        memory: "~{memory_gb} GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/picard:3.1.1--hdfd78af_0"
        maxRetries: 1
   }
}

task quality_score_distribution {
    meta {
        description: "Runs `picard QualityScoreDistribution` to calculate the range of quality scores and creates an accompanying chart"
        external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360037057312-QualityScoreDistribution-Picard-"
        outputs: {
            quality_score_distribution_txt: "The text file output of `QualityScoreDistribution`",
            quality_score_distribution_pdf: "The PDF file output of `QualityScoreDistribution`",
       }
   }

    parameter_meta {
        bam: "Input BAM format file for which to calculate quality score distribution"
        prefix: "Prefix for the output report files. The extensions `.txt` and `.pdf` will be added."
        validation_stringency: {
            description: "Validation stringency for parsing the input BAM.",
            choices: [
                "STRICT",
                "LENIENT",
                "SILENT"
            ],
            tool_default: "STRICT",
       }
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
   }

    input {
        File bam
        String prefix = basename(bam, ".bam") + ".QualityScoreDistribution"
        String validation_stringency = "SILENT"
        Int memory_gb = 8
        Int modify_disk_size_gb = 0
   }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size) + 10 + modify_disk_size_gb
    Int java_heap_size = ceil(memory_gb * 0.9)

    command <<<
        picard -Xmx~{java_heap_size}g QualityScoreDistribution \
            --VALIDATION_STRINGENCY ~{validation_stringency} \
            -I ~{bam} \
            -O ~{prefix}.txt \
            -CHART ~{prefix}.pdf
    >>>

    output {
        File quality_score_distribution_txt = prefix + ".txt"
        File quality_score_distribution_pdf = prefix + ".pdf"
   }

    runtime {
        memory: "~{memory_gb} GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/picard:3.1.1--hdfd78af_0"
        maxRetries: 1
   }
}