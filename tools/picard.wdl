## [Homepage](https://broadinstitute.github.io/picard/)
## 
## This file has been split into two separate files for better organization:
## - picard-qc.wdl: Contains tasks related to QC reporting
## - picard-manipulation.wdl: Contains tasks related to BAM manipulation

version 1.1

import "picard-qc.wdl" as picard_qc
import "picard-manipulation.wdl" as picard_manipulation

# This workflow serves as a wrapper to provide backward compatibility
# with existing pipelines that import tasks from picard.wdl
workflow picard {
    # QC Tasks
    
    # validate_bam task wrapper
    call picard_qc.validate_bam
    
    # collect_wgs_metrics task wrapper
    call picard_qc.collect_wgs_metrics
    
    # collect_alignment_summary_metrics task wrapper
    call picard_qc.collect_alignment_summary_metrics
    
    # collect_gc_bias_metrics task wrapper
    call picard_qc.collect_gc_bias_metrics
    
    # collect_insert_size_metrics task wrapper
    call picard_qc.collect_insert_size_metrics
    
    # quality_score_distribution task wrapper
    call picard_qc.quality_score_distribution
    
    # Manipulation Tasks
    
    # mark_duplicates task wrapper
    call picard_manipulation.mark_duplicates
    
    # sort task wrapper
    call picard_manipulation.sort
    
    # merge_sam_files task wrapper
    call picard_manipulation.merge_sam_files
    
    # clean_sam task wrapper
    call picard_manipulation.clean_sam
    
    # bam_to_fastq task wrapper
    call picard_manipulation.bam_to_fastq
    
    # merge_vcfs task wrapper
    call picard_manipulation.merge_vcfs
    
    # scatter_interval_list task wrapper
    call picard_manipulation.scatter_interval_list
    
    # create_sequence_dictionary task wrapper
    call picard_manipulation.create_sequence_dictionary
    
    output {
        # QC task outputs
        File? validate_report = validate_bam.validate_report
        
        File? wgs_metrics = collect_wgs_metrics.wgs_metrics
        
        File? alignment_metrics = collect_alignment_summary_metrics.alignment_metrics
        File? alignment_metrics_pdf = collect_alignment_summary_metrics.alignment_metrics_pdf
        
        File? gc_bias_metrics = collect_gc_bias_metrics.gc_bias_metrics
        File? gc_bias_metrics_summary = collect_gc_bias_metrics.gc_bias_metrics_summary
        File? gc_bias_metrics_pdf = collect_gc_bias_metrics.gc_bias_metrics_pdf
        
        File? insert_size_metrics = collect_insert_size_metrics.insert_size_metrics
        File? insert_size_metrics_pdf = collect_insert_size_metrics.insert_size_metrics_pdf
        
        File? quality_score_distribution_txt = quality_score_distribution.quality_score_distribution_txt
        File? quality_score_distribution_pdf = quality_score_distribution.quality_score_distribution_pdf
        
        # Manipulation task outputs
        File? duplicate_marked_bam = mark_duplicates.duplicate_marked_bam
        File? duplicate_marked_bam_index = mark_duplicates.duplicate_marked_bam_index
        File? duplicate_marked_bam_md5 = mark_duplicates.duplicate_marked_bam_md5
        File? mark_duplicates_metrics = mark_duplicates.mark_duplicates_metrics
        
        File? sorted_bam = sort.sorted_bam
        File? sorted_bam_index = sort.sorted_bam_index
        File? sorted_bam_md5 = sort.sorted_bam_md5
        
        File? merged_bam = merge_sam_files.merged_bam
        File? merged_bam_index = merge_sam_files.merged_bam_index
        File? merged_bam_md5 = merge_sam_files.merged_bam_md5
        
        File? cleaned_bam = clean_sam.cleaned_bam
        File? cleaned_bam_index = clean_sam.cleaned_bam_index
        File? cleaned_bam_md5 = clean_sam.cleaned_bam_md5
        
        File? read_one_fastq_gz = bam_to_fastq.read_one_fastq_gz
        File? read_two_fastq_gz = bam_to_fastq.read_two_fastq_gz
        
        File? merged_vcf = merge_vcfs.merged_vcf
        File? merged_vcf_index = merge_vcfs.merged_vcf_index
        
        Array[File]? interval_lists_scatter = scatter_interval_list.interval_lists_scatter
        Int? interval_count = scatter_interval_list.interval_count
        
        File? dictionary = create_sequence_dictionary.dictionary
    }
}
