## [Homepage](https://broadinstitute.github.io/picard/)
# SPDX-License-Identifier: MIT

version 1.1

task mark_duplicates {
    meta {
        description: "This WDL task marks duplicate reads in the input BAM file using Picard."
        external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-"
        outputs: {
            duplicate_marked_bam: "The input BAM with computationally determined duplicates marked."
            duplicate_marked_bam_index: "The `.bai` BAM index file associated with `duplicate_marked_bam`"
            duplicate_marked_bam_md5: "The md5sum of `duplicate_marked_bam`"
            mark_duplicates_metrics: {
                description: "The METRICS_FILE result of `picard MarkDuplicates`"
                external_help: "http://broadinstitute.github.io/picard/picard-metric-definitions.html#DuplicationMetrics"
            }
        }
    }

    parameter_meta {
        bam: "Input BAM format file in which to mark duplicates"
        prefix: "Prefix for the MarkDuplicates result files. The extensions `.bam`, `.bam.bai`, `.bam.md5`, and `.metrics.txt` will be added."
        create_bam: "Enable BAM creation (true)? Or only output MarkDuplicates metrics (false)?"
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
        max_retries: "Number of times to retry in case of failure"
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

        picard -Xmx~{java_heap_size}g MarkDuplicates \
            -I ~{bam} \
            -O ~{if create_bam then prefix + ".bam" else "/dev/null"} \
            --VALIDATION_STRINGENCY SILENT \
            --CREATE_INDEX ~{create_bam} \
            --CREATE_MD5_FILE ~{create_bam} \
            --COMPRESSION_LEVEL 5 \
            --METRICS_FILE ~{prefix}.metrics.txt
        
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
        memory: "~{memory_gb} GB"
        disk: "~{disk_size_gb} GB"
        docker: 'quay.io/biocontainers/picard:3.1.0--hdfd78af_0'
        maxRetries: max_retries
    }
}

task validate_bam {
    # TODO should this be refactored to behave as "default" Picard behaves?
    #   Default Picard has some weird/not ideal behaviors
    #   e.g. `max_errors = 100`
    meta {
        description: "This WDL task validates the input BAM file for correct formatting using Picard."
        external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360057440611-ValidateSamFile-Picard-"
        outputs: {
            validate_report: "Validation report produced by `picard ValidateSamFile`. Validation warnings and errors are logged."
            validated_bam: "The unmodified input BAM after it has been succesfully validated"
        }
    }

    parameter_meta {
        bam: "Input BAM format file to validate"
        reference_fasta: "Reference genome in FASTA format. Presence of the reference FASTA allows for `NM` tag validation."
        ignore_list: {
            description: "List of Picard errors and warnings to ignore. Possible values can be found on the GATK website (see `external_help`)."
            external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360035891231-Errors-in-SAM-or-BAM-files-can-be-diagnosed-with-ValidateSamFile"
        }
        outfile_name: "Name for the ValidateSamFile report file"
        succeed_on_errors: "Succeed the task even if errors *and/or* warnings are detected"
        succeed_on_warnings: "Succeed the task if warnings are detected and there are no errors. Overridden by `succeed_on_errors`"
        summary_mode: "Enable SUMMARY mode?"
        index_validation_stringency_less_exhaustive: "Set `INDEX_VALIDATION_STRINGENCY=LESS_EXHAUSTIVE`?"
        max_errors: "Set the value of MAX_OUTPUT for `picard ValidateSamFile`. The Picard default is 100, a lower number can enable fast fail behavior"
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
        max_retries: "Number of times to retry in case of failure"
    }

    input {
        File bam
        File? reference_fasta
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

    String reference_arg = if defined(reference_fasta)
        then "-R ~{reference_fasta}"
        else ""
    String mode_arg = if (summary_mode) then "--MODE SUMMARY" else ""
    String stringency_arg = if (index_validation_stringency_less_exhaustive)
        then "--INDEX_VALIDATION_STRINGENCY LESS_EXHAUSTIVE"
        else ""

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
            echo "Problems detected by Picard ValidateSamFile" > /dev/stderr
            grep -E "$GREP_PATTERN" ~{outfile_name} > /dev/stderr
            exit $rc
        fi
    >>>

    output {
        File validate_report = outfile_name
        File validated_bam = bam
    }

    runtime {
        memory: "~{memory_gb} GB"
        disk: "~{disk_size_gb} GB"
        docker: 'quay.io/biocontainers/picard:3.1.0--hdfd78af_0'
        maxRetries: max_retries
    }
}

task sort {
    meta {
        description: "This WDL task sorts the input BAM file."
        external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360036510732-SortSam-Picard-"
        outputs: {
            sorted_bam: "The input BAM after it has been sorted according to `sort_order`"
            sorted_bam_index: "The `.bai` BAM index file associated with `sorted_bam`"
            sorted_bam_md5: "The md5sum of `sorted_bam`"
        }
    }

    parameter_meta {
        bam: "Input BAM format file to sort"
        sort_order: {
            description: "Order by which to sort the input BAM"
            choices: [
                'queryname',
                'coordinate',
                'duplicate'
            ]
        }
        prefix: "Prefix for the sorted BAM file and accessory files. The extensions `.bam`, `.bam.bai`, and `.bam.md5` will be added."
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
        max_retries: "Number of times to retry in case of failure"
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
        set -euo pipefail

        picard -Xmx~{java_heap_size}g SortSam \
            -I ~{bam} \
            -O ~{outfile_name} \
            -SO ~{sort_order} \
            --CREATE_INDEX true \
            --CREATE_MD5_FILE true \
            --COMPRESSION_LEVEL 5 \
            --VALIDATION_STRINGENCY SILENT
        
        mv ~{prefix}.bai ~{outfile_name}.bai
    >>>

    output {
        File sorted_bam = outfile_name
        File sorted_bam_index = outfile_name + ".bai"
        File sorted_bam_md5 = outfile_name + ".md5"
    }

    runtime {
        memory: "~{memory_gb} GB"
        disk: "~{disk_size_gb} GB"
        docker: 'quay.io/biocontainers/picard:3.1.0--hdfd78af_0'
        maxRetries: max_retries
    }
}

task merge_sam_files {
    meta {
        description: "This WDL task merges the input BAM files into a single BAM file. All input BAMs are assumed to be sorted according to `sort_order`."
        external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360057440751-MergeSamFiles-Picard-"
        outputs: {
            merged_bam: "The BAM resulting from merging all the input BAMs"
            merged_bam_index: "The `.bai` BAM index file associated with `merged_bam`"
            merged_bam_md5: "The md5sum of `merged_bam`"
        }
    }

    parameter_meta {
        bams: "Input BAMs to merge. All BAMs are assumed to be sorted according to `sort_order`."
        prefix: "Prefix for the merged BAM file and accessory files. The extensions `.bam`, `.bam.bai`, and `.bam.md5` will be added."
        sort_order: {
            description: "Sort order for the output merged BAM. It is assumed all input BAMs share this order."
            choices: [
                'unsorted',
                'queryname',
                'coordinate',
                'duplicate',
                'unknown'  # TODO what does this mean?
            ]
        }
        threading: "Option to create a background thread to encode, compress and write to disk the output file. The threaded version uses about 20% more CPU and decreases runtime by ~20% when writing out a compressed BAM file."  # TODO do we need >1 ncpu to enable this?
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
        max_retries: "Number of times to retry in case of failure"
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

    Array[String] input_arg = prefix("--INPUT ", bams)

    String outfile_name = prefix + ".bam"

    command <<<
        set -euo pipefail

        picard -Xmx~{java_heap_size}g MergeSamFiles \
            ~{sep(" ", input_arg)} \
            --OUTPUT ~{outfile_name} \
            --ASSUME_SORTED true \
            --SORT_ORDER ~{sort_order} \
            --USE_THREADING ~{threading} \
            --CREATE_INDEX true \
            --CREATE_MD5_FILE true \
            --VALIDATION_STRINGENCY SILENT
        
        mv ~{prefix}.bai ~{outfile_name}.bai
    >>>

    runtime{
        memory: "~{memory_gb} GB"
        disk: "~{disk_size_gb} GB"
        docker: 'quay.io/biocontainers/picard:3.1.0--hdfd78af_0'
        maxRetries: max_retries
    }

    output {
        File merged_bam = outfile_name
        File merged_bam_index = outfile_name + ".bai"
        File merged_bam_md5 = outfile_name + ".md5"
    }
}

task clean_sam {
    meta {
        description: "This WDL task cleans the input BAM file. Cleans soft-clipping beyond end-of-reference, sets MAPQ=0 for unmapped reads."
        external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360036885571-CleanSam-Picard-"
        outputs: {
            cleaned_bam: "A cleaned version of the input BAM"
            cleaned_bam_index: "The `.bai` BAM index file associated with `cleaned_bam`"
            cleaned_bam_md5: "The md5sum of `cleaned_bam`"
        }
    }

    parameter_meta {
        bam: "Input BAM format file to clean"
        prefix: "Prefix for the cleaned BAM file and accessory files. The extensions `.bam`, `.bam.bai`, and `.bam.md5` will be added."
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
        max_retries: "Number of times to retry in case of failure"
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
        set -euo pipefail

        picard -Xmx~{java_heap_size}g CleanSam \
            -I ~{bam} \
            --CREATE_INDEX true \
            --CREATE_MD5_FILE true \
            -O ~{outfile_name}
        
        mv ~{prefix}.bai ~{outfile_name}.bai
    >>>

    output {
        File cleaned_bam = outfile_name
        File cleaned_bam_index = outfile_name + ".bai"
        File cleaned_bam_md5 = outfile_name + ".md5"
    }

    runtime {
        memory: "~{memory_gb} GB"
        disk: "~{disk_size_gb} GB"
        docker: 'quay.io/biocontainers/picard:3.1.0--hdfd78af_0'
        maxRetries: max_retries
    }
}

task collect_wgs_metrics {
    # TODO not all options exposed
    meta {
        description: "This WDL task runs `picard CollectWgsMetrics`  to collect metrics about the fractions of reads that pass base- and mapping-quality filters as well as coverage (read-depth) levels."
        external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360037226132-CollectWgsMetrics-Picard-"
        outputs: {
            wgs_metrics: {
                description: "Output report of `picard CollectWgsMetrics`"
                external_help: "https://broadinstitute.github.io/picard/picard-metric-definitions.html#CollectWgsMetrics.WgsMetrics"
            }
        }
    }

    parameter_meta {
        bam: "Input BAM format file for which to calculate WGS metrics"
        reference_fasta: "Gzipped reference genome in FASTA format"
        outfile_name: "Name for the metrics result file"
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
        max_retries: "Number of times to retry in case of failure"
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
        memory: "~{memory_gb} GB"
        disk: "~{disk_size_gb} GB"
        docker: 'quay.io/biocontainers/picard:3.1.0--hdfd78af_0'
        maxRetries: max_retries
    }
}

task collect_alignment_summary_metrics {
    # TODO check for other options
    meta {
        description: "This WDL task runs `picard CollectAlignmentSummaryMetrics` to calculate metrics detailing the quality of the read alignments as well as the proportion of the reads that passed machine signal-to-noise threshold quality filters."
        external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360040507751-CollectAlignmentSummaryMetrics-Picard-"
        outputs: {
            alignment_metrics: {
                description: "The text file output of `CollectAlignmentSummaryMetrics`"
                external_help: "http://broadinstitute.github.io/picard/picard-metric-definitions.html#AlignmentSummaryMetrics"
            }
            alignment_metrics_pdf: "The PDF file output of `CollectAlignmentSummaryMetrics`"
        }
    }

    parameter_meta {
        bam: "Input BAM format file for which to calculate alignment metrics"
        prefix: "Prefix for the output report files. The extensions `.txt` and `.pdf` will be added."
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
        max_retries: "Number of times to retry in case of failure"
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
        memory: "~{memory_gb} GB"
        disk: "~{disk_size_gb} GB"
        docker: 'quay.io/biocontainers/picard:3.1.0--hdfd78af_0'
        maxRetries: max_retries
    }
}

task collect_gc_bias_metrics {
    # TODO check for other options
    meta {
        description: "This WDL task runs `picard CollectGcBiasMetrics` to collect information about the relative proportions of guanine (G) and cytosine (C) nucleotides."
        external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360037593931-CollectGcBiasMetrics-Picard-"
        outputs: {
            gc_bias_metrics: {
                description: "The full text file output of `CollectGcBiasMetrics`"
                external_help: "http://broadinstitute.github.io/picard/picard-metric-definitions.html#GcBiasDetailMetrics"
            }
            gc_bias_metrics_summary: {
                description: "The summary text file output of `CollectGcBiasMetrics`"
                external_help: "http://broadinstitute.github.io/picard/picard-metric-definitions.html#GcBiasSummaryMetrics"
            }
            gc_bias_metrics_pdf: "The PDF file output of `CollectGcBiasMetrics`"
        }
    }

    parameter_meta {
        bam: "Input BAM format file for which to calculate GC bias metrics"
        reference_fasta: "Reference sequences in FASTA format"
        prefix: "Prefix for the output report files. The extensions `.txt`, `.summary.txt`, and `.pdf` will be added."
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
        max_retries: "Number of times to retry in case of failure"
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
        memory: "~{memory_gb} GB"
        disk: "~{disk_size_gb} GB"
        docker: 'quay.io/biocontainers/picard:3.1.0--hdfd78af_0'
        maxRetries: max_retries
    }
}

task collect_insert_size_metrics {
    # TODO check for other options
    # TODO what happens if a SE BAM is supplied?
    meta {
        description: "This WDL task runs `picard CollectInsertSizeMetrics` to collect metrics for validating library construction including the insert size distribution and read orientation of paired-end libraries."
        external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360037055772-CollectInsertSizeMetrics-Picard-"
        outputs: {
            insert_size_metrics: {
                description: "The text file output of `CollectInsertSizeMetrics`"
                external_help: "http://broadinstitute.github.io/picard/picard-metric-definitions.html#InsertSizeMetrics"
            }
            insert_size_metrics_pdf: "The PDF file output of `CollectInsertSizeMetrics`"
        }
    }

    parameter_meta {
        bam: "Input BAM format file for which to calculate insert size metrics"
        prefix: "Prefix for the output report files. The extensions `.txt` and `.pdf` will be added."
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
        max_retries: "Number of times to retry in case of failure"
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
        memory: "~{memory_gb} GB"
        disk: "~{disk_size_gb} GB"
        docker: 'quay.io/biocontainers/picard:3.1.0--hdfd78af_0'
        maxRetries: max_retries
    }
}

task quality_score_distribution {
    # TODO check for other options
    meta {
        description: "This WDL task runs `picard QualityScoreDistribution` to calculate the range of quality scores and creates an accompanying chart."
        external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360037057312-QualityScoreDistribution-Picard-"
        outputs: {
            quality_score_distribution_txt: "The text file output of `QualityScoreDistribution`"
            quality_score_distribution_pdf: "The PDF file output of `QualityScoreDistribution`"
        }
    }

    parameter_meta {
        bam: "Input BAM format file for which to calculate quality score distribution"
        prefix: "Prefix for the output report files. The extensions `.txt` and `.pdf` will be added."
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
        max_retries: "Number of times to retry in case of failure"
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
        memory: "~{memory_gb} GB"
        disk: "~{disk_size_gb} GB"
        docker: 'quay.io/biocontainers/picard:3.1.0--hdfd78af_0'
        maxRetries: max_retries
    }
}

task bam_to_fastq {
    meta {
        description: "**[Deprecated]** This WDL task converts the input BAM file into FASTQ format files. This task has been deprecated in favor of `samtools.collate_to_fastq` which is more performant and doesn't error on 'illegal mate states'."
        deprecated: true
    }

    parameter_meta {
        bam: "Input BAM format file to convert to FASTQ"
        paired: "Is the data paired-end (true) or single-end (false)?"
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
        max_retries: "Number of times to retry in case of failure"
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
        memory: "~{memory_gb} GB"
        disk: "~{disk_size_gb} GB"
        docker: 'quay.io/biocontainers/picard:2.27.5--hdfd78af_0'
        maxRetries: max_retries
    }
}
