## [Homepage](https://broadinstitute.github.io/picard/)
## 
## This file has been split into two separate files for better organization:
## - picard-qc.wdl: Contains tasks related to QC reporting
## - picard-manipulation.wdl: Contains tasks related to BAM manipulation

version 1.1

import "picard-qc.wdl" as picard_qc
import "picard-manipulation.wdl" as picard_manipulation

# QC Tasks - wrappers to maintain backward compatibility
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
    Int max_errors = 2147483647
    Int memory_gb = 16
    Int modify_disk_size_gb = 0
  }

  call picard_qc.validate_bam as wrapped {
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
    File validate_report = wrapped.validate_report
  }

  runtime {
    container: wrapped.runtime.container
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

  call picard_qc.collect_wgs_metrics as wrapped {
    input:
      bam = bam,
      reference_fasta = reference_fasta,
      outfile_name = outfile_name,
      validation_stringency = validation_stringency,
      memory_gb = memory_gb,
      modify_disk_size_gb = modify_disk_size_gb
  }

  output {
    File wgs_metrics = wrapped.wgs_metrics
  }

  runtime {
    container: wrapped.runtime.container
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

  call picard_qc.collect_alignment_summary_metrics as wrapped {
    input:
      bam = bam,
      prefix = prefix,
      validation_stringency = validation_stringency,
      memory_gb = memory_gb,
      modify_disk_size_gb = modify_disk_size_gb
  }

  output {
    File alignment_metrics = wrapped.alignment_metrics
    File alignment_metrics_pdf = wrapped.alignment_metrics_pdf
  }

  runtime {
    container: wrapped.runtime.container
  }
}

task collect_gc_bias_metrics {
  input {
    File bam
    File reference_fasta
    String prefix = basename(bam, ".bam") + ".CollectGcBiasMetrics"
    String validation_stringency = "SILENT"
    Int memory_gb = 16
    Int modify_disk_size_gb = 0
  }

  call picard_qc.collect_gc_bias_metrics as wrapped {
    input:
      bam = bam,
      reference_fasta = reference_fasta,
      prefix = prefix,
      validation_stringency = validation_stringency,
      memory_gb = memory_gb,
      modify_disk_size_gb = modify_disk_size_gb
  }

  output {
    File gc_bias_metrics = wrapped.gc_bias_metrics
    File gc_bias_metrics_chart = wrapped.gc_bias_metrics_chart
    File gc_bias_metrics_summary = wrapped.gc_bias_metrics_summary
  }

  runtime {
    container: wrapped.runtime.container
  }
}

task collect_insert_size_metrics {
  input {
    File bam
    String prefix = basename(bam, ".bam") + ".CollectInsertSizeMetrics"
    String validation_stringency = "SILENT"
    Boolean include_duplicates = false
    Int memory_gb = 16
    Int modify_disk_size_gb = 0
  }

  call picard_qc.collect_insert_size_metrics as wrapped {
    input:
      bam = bam,
      prefix = prefix,
      validation_stringency = validation_stringency,
      include_duplicates = include_duplicates,
      memory_gb = memory_gb,
      modify_disk_size_gb = modify_disk_size_gb
  }

  output {
    File insert_size_metrics = wrapped.insert_size_metrics
    File insert_size_histogram = wrapped.insert_size_histogram
  }

  runtime {
    container: wrapped.runtime.container
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

  call picard_qc.quality_score_distribution as wrapped {
    input:
      bam = bam,
      prefix = prefix,
      validation_stringency = validation_stringency,
      memory_gb = memory_gb,
      modify_disk_size_gb = modify_disk_size_gb
  }

  output {
    File quality_score_metrics = wrapped.quality_score_metrics
    File quality_score_chart = wrapped.quality_score_chart
  }

  runtime {
    container: wrapped.runtime.container
  }
}

# Manipulation Tasks - wrappers to maintain backward compatibility
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

  call picard_manipulation.mark_duplicates as wrapped {
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
    File? duplicate_marked_bam = wrapped.duplicate_marked_bam
    File? duplicate_marked_bam_index = wrapped.duplicate_marked_bam_index
    File? duplicate_marked_bam_md5 = wrapped.duplicate_marked_bam_md5
    File mark_duplicates_metrics = wrapped.mark_duplicates_metrics
  }

  runtime {
    container: wrapped.runtime.container
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

  call picard_manipulation.sort as wrapped {
    input:
      bam = bam,
      sort_order = sort_order,
      prefix = prefix,
      validation_stringency = validation_stringency,
      memory_gb = memory_gb,
      modify_disk_size_gb = modify_disk_size_gb
  }

  output {
    File sorted_bam = wrapped.sorted_bam
    File? sorted_bam_index = wrapped.sorted_bam_index
    File sorted_bam_md5 = wrapped.sorted_bam_md5
  }

  runtime {
    container: wrapped.runtime.container
  }
}

task merge_sam_files {
  input {
    Array[File] bams
    String prefix
    String sort_order = "coordinate" 
    String validation_stringency = "SILENT"
    Boolean threading = true
    Int memory_gb = 8
    Int modify_disk_size_gb = 0
  }

  call picard_manipulation.merge_sam_files as wrapped {
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
    File merged_bam = wrapped.merged_bam
    File merged_bam_index = wrapped.merged_bam_index
    File merged_bam_md5 = wrapped.merged_bam_md5
  }

  runtime {
    container: wrapped.runtime.container
  }
}

task clean_sam {
  input {
    File bam
    String prefix = basename(bam, ".bam") + ".CleanSam"
    String validation_stringency = "SILENT"
    Int memory_gb = 8
    Int modify_disk_size_gb = 0
  }

  call picard_manipulation.clean_sam as wrapped {
    input:
      bam = bam,
      prefix = prefix,
      validation_stringency = validation_stringency,
      memory_gb = memory_gb,
      modify_disk_size_gb = modify_disk_size_gb
  }

  output {
    File cleaned_bam = wrapped.cleaned_bam
  }

  runtime {
    container: wrapped.runtime.container
  }
}

task bam_to_fastq {
  input {
    File bam
    String prefix = basename(bam, ".bam") + ".bamToFastq"
    Boolean paired_end = true
    String validation_stringency = "SILENT"
    Int memory_gb = 16
    Int modify_disk_size_gb = 0
  }

  call picard_manipulation.bam_to_fastq as wrapped {
    input:
      bam = bam,
      prefix = prefix,
      paired_end = paired_end,
      validation_stringency = validation_stringency,
      memory_gb = memory_gb,
      modify_disk_size_gb = modify_disk_size_gb
  }

  output {
    File read1_fastq = wrapped.read1_fastq
    File? read2_fastq = wrapped.read2_fastq
  }

  runtime {
    container: wrapped.runtime.container
  }
}

task merge_vcfs {
  input {
    Array[File] vcfs
    String prefix
    String validation_stringency = "SILENT"
    Int memory_gb = 8
    Int modify_disk_size_gb = 0
  }

  call picard_manipulation.merge_vcfs as wrapped {
    input:
      vcfs = vcfs,
      prefix = prefix,
      validation_stringency = validation_stringency,
      memory_gb = memory_gb,
      modify_disk_size_gb = modify_disk_size_gb
  }

  output {
    File merged_vcf = wrapped.merged_vcf
    File merged_vcf_index = wrapped.merged_vcf_index
  }

  runtime {
    container: wrapped.runtime.container
  }
}

task scatter_interval_list {
  input {
    File interval_list
    Int scatter_count
    String prefix = "scattered"
    Int memory_gb = 8
    Int modify_disk_size_gb = 0
  }

  call picard_manipulation.scatter_interval_list as wrapped {
    input:
      interval_list = interval_list,
      scatter_count = scatter_count,
      prefix = prefix,
      memory_gb = memory_gb,
      modify_disk_size_gb = modify_disk_size_gb
  }

  output {
    Array[File] scattered_interval_lists = wrapped.scattered_interval_lists
  }

  runtime {
    container: wrapped.runtime.container
  }
}

task create_sequence_dictionary {
  input {
    File reference_fasta
    String? outfile_name
    Int memory_gb = 8
    Int modify_disk_size_gb = 0
  }

  call picard_manipulation.create_sequence_dictionary as wrapped {
    input:
      reference_fasta = reference_fasta,
      outfile_name = outfile_name,
      memory_gb = memory_gb,
      modify_disk_size_gb = modify_disk_size_gb
  }

  output {
    File sequence_dictionary = wrapped.sequence_dictionary
  }

  runtime {
    container: wrapped.runtime.container
  }
}
