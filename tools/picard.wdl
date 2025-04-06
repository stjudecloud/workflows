## [Homepage](https://broadinstitute.github.io/picard/)
## 
## This file has been split into two separate files for better organization:
## - picard-qc.wdl: Contains tasks related to QC reporting
## - picard-manipulation.wdl: Contains tasks related to BAM manipulation

version 1.1

import "picard-qc.wdl" as picard_qc
import "picard-manipulation.wdl" as picard_manipulation

# QC Tasks - aliased at namespace level for backward compatibility
alias picard_qc.validate_bam as validate_bam
alias picard_qc.collect_wgs_metrics as collect_wgs_metrics
alias picard_qc.collect_alignment_summary_metrics as collect_alignment_summary_metrics
alias picard_qc.collect_gc_bias_metrics as collect_gc_bias_metrics
alias picard_qc.collect_insert_size_metrics as collect_insert_size_metrics
alias picard_qc.quality_score_distribution as quality_score_distribution

# Manipulation Tasks - aliased at namespace level for backward compatibility
alias picard_manipulation.mark_duplicates as mark_duplicates
alias picard_manipulation.sort as sort
alias picard_manipulation.merge_sam_files as merge_sam_files
alias picard_manipulation.clean_sam as clean_sam
alias picard_manipulation.bam_to_fastq as bam_to_fastq
alias picard_manipulation.merge_vcfs as merge_vcfs
alias picard_manipulation.scatter_interval_list as scatter_interval_list
alias picard_manipulation.create_sequence_dictionary as create_sequence_dictionary
