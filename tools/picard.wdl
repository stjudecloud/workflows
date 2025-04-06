## [Homepage](https://broadinstitute.github.io/picard/)
## 
## This file has been split into two separate files for better organization:
## - picard-qc.wdl: Contains tasks related to QC reporting
## - picard-manipulation.wdl: Contains tasks related to BAM manipulation

version 1.1

import "picard-qc.wdl" as picard_qc
import "picard-manipulation.wdl" as picard_manipulation

# QC Tasks
task validate_bam {
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

    call picard_qc.validate_bam {
        input:
            bam = bam,
            reference_fasta = reference_fasta,
            ignore_list = ignore_list, 
            outfile_name = outfile_name,
            validation_stringency = validation_stringency,
            succeed_on_errors = succeed_on_errors,
            succeed_on_warnings = succeed_on_warnings,
            summary_mode = summary_mode,
            index_validation_stringency_less_exhaustive = index_validation_stringency_less_exhaustive,
            max_errors = max_errors,
            memory_gb = memory_gb,
            modify_disk_size_gb = modify_disk_size_gb
    }

    output {
        File validate_report = picard_qc.validate_bam.validate_report
    }
}

task collect_wgs_metrics {
    input {
        File bam
        File reference_fasta
        String outfile_name = basename(bam, ".bam") + ".CollectWgsMetrics.txt"
        String validation_stringency = "SILENT"
        Int memory_gb = 12
        Int modify_disk_size_gb = 0
    }

    call picard_qc.collect_wgs_metrics {
        input:
            bam = bam,
            reference_fasta = reference_fasta,
            outfile_name = outfile_name,
            validation_stringency = validation_stringency,
            memory_gb = memory_gb,
            modify_disk_size_gb = modify_disk_size_gb
    }

    output {
        File wgs_metrics = picard_qc.collect_wgs_metrics.wgs_metrics
    }
}

task collect_alignment_summary_metrics {
    input {
        File bam
        String prefix = basename(bam, ".bam") + ".CollectAlignmentSummaryMetrics"
        String validation_stringency = "SILENT"
        Int memory_gb = 8
        Int modify_disk_size_gb = 0
    }

    call picard_qc.collect_alignment_summary_metrics {
        input:
            bam = bam,
            prefix = prefix,
            validation_stringency = validation_stringency,
            memory_gb = memory_gb,
            modify_disk_size_gb = modify_disk_size_gb
    }

    output {
        File alignment_metrics = picard_qc.collect_alignment_summary_metrics.alignment_metrics
        File alignment_metrics_pdf = picard_qc.collect_alignment_summary_metrics.alignment_metrics_pdf
    }
}

task collect_gc_bias_metrics {
    input {
        File bam
        File reference_fasta
        String prefix = basename(bam, ".bam") + ".CollectGcBiasMetrics"
        String validation_stringency = "SILENT"
        Int memory_gb = 8
        Int modify_disk_size_gb = 0
    }

    call picard_qc.collect_gc_bias_metrics {
        input:
            bam = bam,
            reference_fasta = reference_fasta,
            prefix = prefix,
            validation_stringency = validation_stringency,
            memory_gb = memory_gb,
            modify_disk_size_gb = modify_disk_size_gb
    }

    output {
        File gc_bias_metrics = picard_qc.collect_gc_bias_metrics.gc_bias_metrics
        File gc_bias_metrics_summary = picard_qc.collect_gc_bias_metrics.gc_bias_metrics_summary
        File gc_bias_metrics_pdf = picard_qc.collect_gc_bias_metrics.gc_bias_metrics_pdf
    }
}

task collect_insert_size_metrics {
    input {
        File bam
        String prefix = basename(bam, ".bam") + ".CollectInsertSizeMetrics"
        String validation_stringency = "SILENT"
        Int memory_gb = 8
        Int modify_disk_size_gb = 0
    }

    call picard_qc.collect_insert_size_metrics {
        input:
            bam = bam,
            prefix = prefix,
            validation_stringency = validation_stringency,
            memory_gb = memory_gb,
            modify_disk_size_gb = modify_disk_size_gb
    }

    output {
        File insert_size_metrics = picard_qc.collect_insert_size_metrics.insert_size_metrics
        File insert_size_metrics_pdf = picard_qc.collect_insert_size_metrics.insert_size_metrics_pdf
    }
}

task quality_score_distribution {
    input {
        File bam
        String prefix = basename(bam, ".bam") + ".QualityScoreDistribution"
        String validation_stringency = "SILENT"
        Int memory_gb = 8
        Int modify_disk_size_gb = 0
    }

    call picard_qc.quality_score_distribution {
        input:
            bam = bam,
            prefix = prefix,
            validation_stringency = validation_stringency,
            memory_gb = memory_gb,
            modify_disk_size_gb = modify_disk_size_gb
    }

    output {
        File quality_score_distribution_txt = picard_qc.quality_score_distribution.quality_score_distribution_txt
        File quality_score_distribution_pdf = picard_qc.quality_score_distribution.quality_score_distribution_pdf
    }
}

# Manipulation Tasks
task mark_duplicates {
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

    call picard_manipulation.mark_duplicates {
        input:
            bam = bam,
            prefix = prefix,
            duplicate_scoring_strategy = duplicate_scoring_strategy,
            read_name_regex = read_name_regex,
            tagging_policy = tagging_policy,
            validation_stringency = validation_stringency,
            create_bam = create_bam,
            clear_dt = clear_dt,
            remove_duplicates = remove_duplicates,
            remove_sequencing_duplicates = remove_sequencing_duplicates,
            optical_distance = optical_distance,
            modify_memory_gb = modify_memory_gb,
            modify_disk_size_gb = modify_disk_size_gb
    }

    output {
        File? duplicate_marked_bam = picard_manipulation.mark_duplicates.duplicate_marked_bam
        File? duplicate_marked_bam_index = picard_manipulation.mark_duplicates.duplicate_marked_bam_index
        File? duplicate_marked_bam_md5 = picard_manipulation.mark_duplicates.duplicate_marked_bam_md5
        File mark_duplicates_metrics = picard_manipulation.mark_duplicates.mark_duplicates_metrics
    }
}

task sort {
    input {
        File bam
        String sort_order = "coordinate"
        String prefix = basename(bam, ".bam") + ".sorted"
        String validation_stringency = "SILENT"
        Int memory_gb = 25
        Int modify_disk_size_gb = 0
    }

    call picard_manipulation.sort {
        input:
            bam = bam,
            sort_order = sort_order,
            prefix = prefix,
            validation_stringency = validation_stringency,
            memory_gb = memory_gb,
            modify_disk_size_gb = modify_disk_size_gb
    }

    output {
        File sorted_bam = picard_manipulation.sort.sorted_bam
        File? sorted_bam_index = picard_manipulation.sort.sorted_bam_index
        File sorted_bam_md5 = picard_manipulation.sort.sorted_bam_md5
    }
}

task merge_sam_files {
    input {
        Array[File] bams
        String prefix
        String sort_order = "coordinate"
        String validation_stringency = "SILENT"
        Boolean threading = true
        Int memory_gb = 40
        Int modify_disk_size_gb = 0
    }

    call picard_manipulation.merge_sam_files {
        input:
            bams = bams,
            prefix = prefix,
            sort_order = sort_order,
            validation_stringency = validation_stringency,
            threading = threading,
            memory_gb = memory_gb,
            modify_disk_size_gb = modify_disk_size_gb
    }

    output {
        File merged_bam = picard_manipulation.merge_sam_files.merged_bam
        File merged_bam_index = picard_manipulation.merge_sam_files.merged_bam_index
        File merged_bam_md5 = picard_manipulation.merge_sam_files.merged_bam_md5
    }
}

task clean_sam {
    input {
        File bam
        String prefix = basename(bam, ".bam") + ".cleaned"
        String validation_stringency = "SILENT"
        Int memory_gb = 25
        Int modify_disk_size_gb = 0
    }

    call picard_manipulation.clean_sam {
        input:
            bam = bam,
            prefix = prefix,
            validation_stringency = validation_stringency,
            memory_gb = memory_gb,
            modify_disk_size_gb = modify_disk_size_gb
    }

    output {
        File cleaned_bam = picard_manipulation.clean_sam.cleaned_bam
        File cleaned_bam_index = picard_manipulation.clean_sam.cleaned_bam_index
        File cleaned_bam_md5 = picard_manipulation.clean_sam.cleaned_bam_md5
    }
}

task bam_to_fastq {
    input {
        File bam
        String prefix = basename(bam, ".bam")
        Boolean paired = true
        Int memory_gb = 56
        Int modify_disk_size_gb = 0
    }

    call picard_manipulation.bam_to_fastq {
        input:
            bam = bam,
            prefix = prefix,
            paired = paired,
            memory_gb = memory_gb,
            modify_disk_size_gb = modify_disk_size_gb
    }

    output {
        File read_one_fastq_gz = picard_manipulation.bam_to_fastq.read_one_fastq_gz
        File? read_two_fastq_gz = picard_manipulation.bam_to_fastq.read_two_fastq_gz
    }
}

task merge_vcfs {
    input {
        Array[File] vcfs
        Array[File] vcfs_indexes
        String output_vcf_name
        Int modify_disk_size_gb = 0
    }

    call picard_manipulation.merge_vcfs {
        input:
            vcfs = vcfs,
            vcfs_indexes = vcfs_indexes,
            output_vcf_name = output_vcf_name,
            modify_disk_size_gb = modify_disk_size_gb
    }

    output {
        File merged_vcf = picard_manipulation.merge_vcfs.merged_vcf
        File merged_vcf_index = picard_manipulation.merge_vcfs.merged_vcf_index
    }
}

task scatter_interval_list {
    input {
        File interval_list
        Int scatter_count
        String subdivision_mode = "BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW"
        Boolean unique = true
        Boolean sort = true
    }

    call picard_manipulation.scatter_interval_list {
        input:
            interval_list = interval_list,
            scatter_count = scatter_count,
            subdivision_mode = subdivision_mode,
            unique = unique,
            sort = sort
    }

    output {
        Array[File] interval_lists_scatter = picard_manipulation.scatter_interval_list.interval_lists_scatter
        Int interval_count = picard_manipulation.scatter_interval_list.interval_count
    }
}

task create_sequence_dictionary {
    input {
        File fasta
        String? assembly_name
        String? fasta_url
        String? species
        String outfile_name = basename(fasta, ".fa") + ".dict"
        Int memory_gb = 16
        Int modify_disk_size_gb = 0
    }

    call picard_manipulation.create_sequence_dictionary {
        input:
            fasta = fasta,
            assembly_name = assembly_name,
            fasta_url = fasta_url,
            species = species,
            outfile_name = outfile_name,
            memory_gb = memory_gb,
            modify_disk_size_gb = modify_disk_size_gb
    }

    output {
        File dictionary = picard_manipulation.create_sequence_dictionary.dictionary
    }
}
