## [Homepage](https://broadinstitute.github.io/picard/)

version 1.1

task mark_duplicates {
    meta {
        description: "Marks duplicate reads in the input BAM file using Picard"
        external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-"
        help: "For non-primary reads and unmapped mates of duplicate reads to be included in duplicate analysis, input BAM must be collated. See `external_help` for more information."
        outputs: {
            duplicate_marked_bam: "The input BAM with computationally determined duplicates marked.",
            duplicate_marked_bam_index: "The `.bai` BAM index file associated with `duplicate_marked_bam`",
            duplicate_marked_bam_md5: "The md5sum of `duplicate_marked_bam`",
            mark_duplicates_metrics: {
                description: "The METRICS_FILE result of `picard MarkDuplicates`",
                external_help: "http://broadinstitute.github.io/picard/picard-metric-definitions.html#DuplicationMetrics",
            },
        }
    }

    parameter_meta {
        bam: "Input BAM format file in which to mark duplicates"
        prefix: "Prefix for the MarkDuplicates result files. The extensions `.bam`, `.bam.bai`, `.bam.md5`, and `.metrics.txt` will be added."
        duplicate_scoring_strategy: {
            description: "Strategy for scoring duplicates.",
            choices: [
                "SUM_OF_BASE_QUALITIES",
                "TOTAL_MAPPED_REFERENCE_LENGTH",
                "RANDOM",
            ],
        }
        read_name_regex: "Regular expression for extracting tile names, x coordinates, and y coordinates from read names. The default works for typical Illumina read names."
        tagging_policy: {
            description: "Tagging policy for the output BAM.",
            choices: [
                "DontTag",
                "OpticalOnly",
                "All",
            ],
        }
        validation_stringency: {
            description: "Validation stringency for parsing the input BAM.",
            choices: [
                "STRICT",
                "LENIENT",
                "SILENT",
            ],
            tool_default: "STRICT",
        }
        create_bam: {
            description: "Enable BAM creation (true)? Or only output MarkDuplicates metrics (false)?",
            group: "Common",
        }
        clear_dt: "Clear the `DT` tag from the input BAM? For increased performance, if the input BAM does not have the `DT` tag, set to `false`."
        remove_duplicates: "Remove duplicate reads from the output BAM? If `true`, the output BAM will not contain any duplicate reads."
        remove_sequencing_duplicates: "Remove sequencing duplicates (i.e. optical duplicates) from the output BAM? If `true`, the output BAM will not contain any sequencing duplicates (optical duplicates)."
        optical_distance: "Maximum distance between read coordinates to consider them optical duplicates. If `0`, then optical duplicate marking is disabled. Suggested settings of 100 for unpatterned versions of the Illumina platform (e.g. HiSeq) or 2500 for patterned flowcell models (e.g. NovaSeq). Calculation of distance depends on coordinate data embedded in the read names, typically produced by the Illumina sequencing machines. Optical duplicate detection will not work on non-standard names without modifying `read_name_regex`."
        modify_memory_gb: "Add to or subtract from the default memory allocation. Default memory allocation is determined by the size of the input BAM. Specified in GB."
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        String prefix = basename(bam, ".bam") + ".MarkDuplicates"
        String duplicate_scoring_strategy = "SUM_OF_BASE_QUALITIES"
        String read_name_regex = "^[!-9;-?A-~:]+:([!-9;-?A-~]+):([0-9]+):([0-9]+)$"
        String tagging_policy = "All"
        String validation_stringency = "SILENT"
        Boolean create_bam = true
        Boolean clear_dt = true
        Boolean remove_duplicates = false
        Boolean remove_sequencing_duplicates = false
        Int optical_distance = 0
        Int modify_memory_gb = 0
        Int modify_disk_size_gb = 0
    }

    Float bam_size = size(bam, "GiB")
    Int memory_gb = min(ceil(bam_size + 6), 50) + modify_memory_gb
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
            -I "~{bam}" \
            --METRICS_FILE "~{prefix}.metrics.txt" \
            -O "~{if create_bam then prefix + ".bam" else "/dev/null"}" \
            --CREATE_INDEX ~{create_bam} \
            --CREATE_MD5_FILE ~{create_bam} \
            --VALIDATION_STRINGENCY "~{validation_stringency}" \
            --DUPLICATE_SCORING_STRATEGY "~{duplicate_scoring_strategy}" \
            --READ_NAME_REGEX '~{
                if (optical_distance > 0) then read_name_regex else "null"
            }' \
            --TAGGING_POLICY "~{tagging_policy}" \
            --CLEAR_DT ~{clear_dt} \
            --REMOVE_DUPLICATES ~{remove_duplicates} \
            --REMOVE_SEQUENCING_DUPLICATES ~{remove_sequencing_duplicates} \
            --OPTICAL_DUPLICATE_PIXEL_DISTANCE ~{optical_distance}

        if ~{create_bam}; then
            mv "~{prefix}.bai" "~{prefix}.bam.bai"
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
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/picard:3.1.1--hdfd78af_0"
        maxRetries: 1
    }
}

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
            group: "Common",
        }
        outfile_name: "Name for the ValidateSamFile report file"
        validation_stringency: {
            description: "Validation stringency for parsing the input BAM.",
            choices: [
                "STRICT",
                "LENIENT",
                "SILENT",
            ],
            tool_default: "STRICT",
        }
        succeed_on_errors: {
            description: "Succeed the task even if errors *and/or* warnings are detected",
            group: "Common",
        }
        succeed_on_warnings: {
            description: "Succeed the task if warnings are detected and there are no errors. Overridden by `succeed_on_errors`",
            group: "Common",
        }
        summary_mode: {
            description: "Enable SUMMARY mode?",
            group: "Common",
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
            -I "~{bam}" \
            ~{reference_arg} \
            ~{mode_arg} \
            ~{stringency_arg} \
            --VALIDATION_STRINGENCY "~{validation_stringency}" \
            ~{sep(" ", prefix("--IGNORE ", squote(ignore_list)))} \
            --MAX_OUTPUT ~{max_errors} \
            > "~{outfile_name}" \
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
            && [ "$(grep -Ec "$GREP_PATTERN" "~{outfile_name}")" -gt 0 ]
        then
            >&2 echo "Problems detected by Picard ValidateSamFile"
            >&2 grep -E "$GREP_PATTERN" "~{outfile_name}"
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

task sort {
    meta {
        description: "Sorts the input BAM file"
        external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360036510732-SortSam-Picard-"
        outputs: {
            sorted_bam: "The input BAM after it has been sorted according to `sort_order`",
            sorted_bam_index: "The `.bai` BAM index file associated with `sorted_bam`",
            sorted_bam_md5: "The md5sum of `sorted_bam`",
        }
    }

    parameter_meta {
        bam: "Input BAM format file to sort"
        sort_order: {
            description: "Order by which to sort the input BAM",
            choices: [
                "queryname",
                "coordinate",
                "duplicate",
            ],
            group: "Common",
        }
        prefix: "Prefix for the sorted BAM file and accessory files. The extensions `.bam`, `.bam.bai`, and `.bam.md5` will be added."
        validation_stringency: {
            description: "Validation stringency for parsing the input BAM.",
            choices: [
                "STRICT",
                "LENIENT",
                "SILENT",
            ],
            tool_default: "STRICT",
        }
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        String sort_order = "coordinate"
        String prefix = basename(bam, ".bam") + ".sorted"
        String validation_stringency = "SILENT"
        Int memory_gb = 25
        Int modify_disk_size_gb = 0
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size * 4) + 10 + modify_disk_size_gb
    Int java_heap_size = ceil(memory_gb * 0.9)

    String outfile_name = prefix + ".bam"

    command <<<
        set -euo pipefail

        picard -Xmx~{java_heap_size}g SortSam \
            -I "~{bam}" \
            -O "~{outfile_name}" \
            -SO "~{sort_order}" \
            --CREATE_INDEX true \
            --CREATE_MD5_FILE true \
            --VALIDATION_STRINGENCY "~{validation_stringency}"

        # CREATE_INDEX true only applies when the sort order
        # is coordinate. So the index may not exist.
        if [ -f "~{prefix}.bai" ]; then
            mv "~{prefix}.bai" "~{outfile_name}.bai"
        fi
    >>>

    output {
        File sorted_bam = outfile_name
        File? sorted_bam_index = outfile_name + ".bai"
        File sorted_bam_md5 = outfile_name + ".md5"
    }

    runtime {
        memory: "~{memory_gb} GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/picard:3.1.1--hdfd78af_0"
        maxRetries: 1
    }
}

task merge_sam_files {
    meta {
        description: "Merges the input BAM files into a single BAM file. All input BAMs are assumed to be sorted according to `sort_order`."
        external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360057440751-MergeSamFiles-Picard-"
        outputs: {
            merged_bam: "The BAM resulting from merging all the input BAMs",
            merged_bam_index: "The `.bai` BAM index file associated with `merged_bam`",
            merged_bam_md5: "The md5sum of `merged_bam`",
        }
    }

    parameter_meta {
        bams: "Input BAMs to merge. All BAMs are assumed to be sorted according to `sort_order`."
        prefix: "Prefix for the merged BAM file and accessory files. The extensions `.bam`, `.bam.bai`, and `.bam.md5` will be added."
        sort_order: {
            description: "Sort order for the output merged BAM. It is assumed all input BAMs share this order.",
            choices: [
                "unsorted",
                "queryname",
                "coordinate",
                "duplicate",
            ],
            group: "Common",
        }
        validation_stringency: {
            description: "Validation stringency for parsing the input BAM.",
            choices: [
                "STRICT",
                "LENIENT",
                "SILENT",
            ],
            tool_default: "STRICT",
        }
        threading: "Option to create a background thread to encode, compress and write to disk the output file. The threaded version uses about 20% more CPU and decreases runtime by ~20% when writing out a compressed BAM file. **Sets `runtime.cpu = 2` if `true`. `runtime.cpu = 1` if `false`.**"
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        Array[File] bams
        String prefix
        String sort_order = "coordinate"
        String validation_stringency = "SILENT"
        Boolean threading = true
        Int memory_gb = 40
        Int modify_disk_size_gb = 0
    }

    Float bams_size = size(bams, "GiB")
    Int disk_size_gb = ceil(bams_size * 2) + 10 + modify_disk_size_gb
    Int java_heap_size = ceil(memory_gb * 0.9)

    Array[String] input_arg = prefix("--INPUT ", squote(bams))

    String outfile_name = prefix + ".bam"

    command <<<
        set -euo pipefail

        picard -Xmx~{java_heap_size}g MergeSamFiles \
            ~{sep(" ", input_arg)} \
            --OUTPUT "~{outfile_name}" \
            --ASSUME_SORTED true \
            --SORT_ORDER "~{sort_order}" \
            --USE_THREADING ~{threading} \
            --CREATE_INDEX true \
            --CREATE_MD5_FILE true \
            --VALIDATION_STRINGENCY "~{validation_stringency}"

        mv "~{prefix}.bai" "~{outfile_name}.bai"
    >>>

    output {
        File merged_bam = outfile_name
        File merged_bam_index = outfile_name + ".bai"
        File merged_bam_md5 = outfile_name + ".md5"
    }

    runtime{
        cpu: if threading then 2 else 1
        memory: "~{memory_gb} GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/picard:3.1.1--hdfd78af_0"
        maxRetries: 1
    }
}

task clean_sam {
    meta {
        description: "Cleans the input BAM file. Cleans soft-clipping beyond end-of-reference, sets MAPQ=0 for unmapped reads."
        external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360036885571-CleanSam-Picard-"
        outputs: {
            cleaned_bam: "A cleaned version of the input BAM",
            cleaned_bam_index: "The `.bai` BAM index file associated with `cleaned_bam`",
            cleaned_bam_md5: "The md5sum of `cleaned_bam`",
        }
    }

    parameter_meta {
        bam: "Input BAM format file to clean"
        prefix: "Prefix for the cleaned BAM file and accessory files. The extensions `.bam`, `.bam.bai`, and `.bam.md5` will be added."
        validation_stringency: {
            description: "Validation stringency for parsing the input BAM.",
            choices: [
                "STRICT",
                "LENIENT",
                "SILENT",
            ],
            tool_default: "STRICT",
        }
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        String prefix = basename(bam, ".bam") + ".cleaned"
        String validation_stringency = "SILENT"
        Int memory_gb = 25
        Int modify_disk_size_gb = 0
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size * 2) + 10 + modify_disk_size_gb
    Int java_heap_size = ceil(memory_gb * 0.9)

    String outfile_name = prefix + ".bam"

    command <<<
        set -euo pipefail

        picard -Xmx~{java_heap_size}g CleanSam \
            -I "~{bam}" \
            --CREATE_INDEX true \
            --CREATE_MD5_FILE true \
            --VALIDATION_STRINGENCY "~{validation_stringency}" \
            -O "~{outfile_name}"

        mv "~{prefix}.bai" "~{outfile_name}.bai"
    >>>

    output {
        File cleaned_bam = outfile_name
        File cleaned_bam_index = outfile_name + ".bai"
        File cleaned_bam_md5 = outfile_name + ".md5"
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
                "SILENT",
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
            -I "~{bam}" \
            -R "~{reference_fasta}" \
            -O "~{outfile_name}" \
            --VALIDATION_STRINGENCY "~{validation_stringency}" \
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
                "SILENT",
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
            -I "~{bam}" \
            --VALIDATION_STRINGENCY "~{validation_stringency}" \
            -O "~{prefix}.txt" \
            -H "~{prefix}.pdf"
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
                "SILENT",
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
            -I "~{bam}" \
            -R "~{reference_fasta}" \
            --VALIDATION_STRINGENCY "~{validation_stringency}" \
            -O "~{prefix}.txt" \
            -S "~{prefix}.summary.txt" \
            -CHART "~{prefix}.pdf"
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
                "SILENT",
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
            -I "~{bam}" \
            --VALIDATION_STRINGENCY "~{validation_stringency}" \
            -O "~{prefix}.txt" \
            -H "~{prefix}.pdf"
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
                "SILENT",
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
            --VALIDATION_STRINGENCY "~{validation_stringency}" \
            -I "~{bam}" \
            -O "~{prefix}.txt" \
            -CHART "~{prefix}.pdf"
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

#@ except: MatchingOutputMeta
task bam_to_fastq {
    meta {
        description: "**[Deprecated]** This WDL task converts the input BAM file into FASTQ format files. This task has been deprecated in favor of `samtools.bam_to_fastq` which is more performant and doesn't error on 'illegal mate states'."
        deprecated: true
    }

    parameter_meta {
        bam: "Input BAM format file to convert to FASTQ"
        prefix: "Prefix for the <type of file> file. The extension `<extension>` will be added."
        paired: {
            description: "Is the data Paired-End (true) or Single-End (false)?",
            group: "Common",
        }
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        String prefix = basename(bam, ".bam")
        Boolean paired = true
        Int memory_gb = 56
        Int modify_disk_size_gb = 0
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size * 4) + 10 + modify_disk_size_gb
    Int java_heap_size = ceil(memory_gb * 0.9)

    command <<<
        set -euo pipefail

        picard -Xmx~{java_heap_size}g SamToFastq INPUT="~{bam}" \
            FASTQ="~{prefix}.R1.fastq" \
            ~{(
                if paired
                then "SECOND_END_FASTQ=" + prefix + ".R2.fastq"
                else ""
            )} \
            RE_REVERSE=true \
            VALIDATION_STRINGENCY=SILENT

        gzip "~{prefix}.R1.fastq" \
            ~{if paired then prefix + ".R2.fastq" else ""}
    >>>

    output {
        File read_one_fastq_gz = "~{prefix}.R1.fastq.gz"
        File? read_two_fastq_gz = "~{prefix}.R2.fastq.gz"
    }

    runtime{
        memory: "~{memory_gb} GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/picard:3.1.1--hdfd78af_0"
        maxRetries: 1
    }
}

task merge_vcfs {
    meta {
        description: "Merges the input VCF files into a single VCF file"
        external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360036713331-MergeVcfs-Picard"
        outputs: {
            merged_vcf: "The merged VCF file",
            merged_vcf_index: "The index file associated with the merged VCF file",
        }
    }

    parameter_meta {
        vcfs: "Input VCF format files to merge. May be gzipped or binary compressed."
        vcfs_indexes: "Index files associated with the input VCF files"
        output_vcf_name: "Name for the merged VCF file"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        Array[File] vcfs
        Array[File] vcfs_indexes
        String output_vcf_name
        Int modify_disk_size_gb = 0
    }

    Int disk_size_gb = ceil(size(vcfs, "GiB") * 2) + 10 + modify_disk_size_gb

    command <<<
        picard -Xms2000m \
            MergeVcfs \
            ~{sep(" ", prefix("--INPUT ", squote(vcfs)))} \
            --OUTPUT "~{output_vcf_name}"
    >>>

    output {
        File merged_vcf = output_vcf_name
        File merged_vcf_index = "~{output_vcf_name}.tbi"
    }

    runtime {
        memory: "4 GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/picard:2.27.5--hdfd78af_0"
        maxRetries: 1
    }
}

task scatter_interval_list {
    meta {
        description: "Splits an interval list into smaller interval lists for parallel processing"
        external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360036897212-IntervalListTools-Picard"
        outputs: {
            interval_lists_scatter: "The split interval lists",
            interval_count: "The number of split interval lists",
        }
    }

    parameter_meta  {
        interval_list: "Input interval list to split"
        scatter_count: "Number of interval lists to create"
        subdivision_mode: {
            description: "How to subdivide the intervals",
            choices: [
                "BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW",
                "INTERVAL_SUBDIVISION",
                "BALANCING_WITHOUT_INTERVAL_SUBDIVISION",
            ],
        }
        unique: "Should the output interval lists contain unique intervals? Implies sort=true. Merges overlapping or adjacent intervals."
        sort: "Should the output interval lists be sorted? Sorts by coordinate."
    }

    input {
        File interval_list
        Int scatter_count
        String subdivision_mode = "BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW"
        Boolean unique = true
        Boolean sort = true
    }

    command <<<
        set -euo pipefail

        mkdir out
        picard -Xms1g \
            IntervalListTools \
            --SCATTER_COUNT ~{scatter_count} \
            --SUBDIVISION_MODE "~{subdivision_mode}" \
            --UNIQUE ~{unique} \
            --SORT ~{sort} \
            --INPUT "~{interval_list}" \
            --OUTPUT out

        bash <<CODE
        I=0
        for list in out/*/*.interval_list
        do
           I=\$((I+1))
           dir=\$(dirname \$list)
           name=\$(basename \$list)
           mv \$list \${dir}/\${I}\${name}
        done
        echo \$I > interval_count.txt
        CODE
    >>>

    output {
        Array[File] interval_lists_scatter = glob("out/*/*.interval_list")
        Int interval_count = read_int("interval_count.txt")
    }

    runtime {
        memory: "2 GB"
        disks: "1 GB"
        container: "quay.io/biocontainers/picard:2.27.5--hdfd78af_0"
        maxRetries: 1
    }
}

task create_sequence_dictionary {
    meta {
        description: "Creates a sequence dictionary for the input FASTA file using Picard"
        external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/13832748622491-CreateSequenceDictionary-Picard-"
        outputs: {
            dictionary: "Sequence dictionary produced by `picard CreateSequenceDictionary`."
        }
    }

    parameter_meta {
        fasta: "Input FASTA format file from which to create dictionary"
        assembly_name: "Value to put in AS field of sequence dictionary"
        fasta_url: "Value to put in UR field of sequence dictionary"
        species: "Value to put in SP field of sequence dictionary"
        outfile_name: "Name for the CreateSequenceDictionary dictionary file"
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File fasta
        String? assembly_name
        String? fasta_url
        String? species
        String outfile_name = basename(fasta, ".fa") + ".dict"
        Int memory_gb = 16
        Int modify_disk_size_gb = 0
    }

    Float fasta_size = size(fasta, "GiB")
    Int disk_size_gb = ceil(fasta_size * 2) + 10 + modify_disk_size_gb
    Int java_heap_size = ceil(memory_gb * 0.9)

    command <<<
        set -euo pipefail

        picard -Xmx~{java_heap_size}g CreateSequenceDictionary \
            -R "~{fasta}" \
            ~{(
                if defined(assembly_name)
                then "--GENOME_ASSEMBLY " + assembly_name
                else ""
            )} \
            ~{if defined(fasta_url) then "--URI " + fasta_url else ""} \
            ~{if defined(species) then "--SPECIES " + species else ""} \
            > "~{outfile_name}"
    >>>

    output {
        File dictionary = outfile_name
    }

    runtime {
        cpu: 1
        memory: "~{memory_gb} GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/picard:3.1.0--hdfd78af_0"
        maxRetries: 1
    }
}
