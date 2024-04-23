# SPDX-License-Identifier: MIT
# Copyright St. Jude Children's Research Hospital
version 1.1

import "../../tools/picard.wdl"
import "../../tools/gatk4.wdl" as gatk

workflow rnaseq_variant_calling {
    meta {
        description: "Call variants from RNA-Seq data. Produces a VCF file of variants. Based on GATK RNA-Seq short variant calling best practices pipeline."
        outputs: {
            recalibrated_bam: "BAM that has undergone recalibration of base quality scores"
            recalibrated_bam_index: "Index file for recalibrated BAM file"
            merged_vcf: "VCF of raw RNA-Seq variant calls"
            merged_vcf_index: "Index for merged VCF file"
            variant_filtered_vcf: "VCF file after variant filters have been applied"
            variant_filtered_vcf_index: "Index for filtered variant VCF file"
        }
    }

    parameter_meta {
        bam: "BAM file of aligned RNA-Seq reads"
        bam_index: "Index file for BAM file"
        fasta: "Reference FASTA file"
        fasta_index: "Index file for reference FASTA file"
        dict: "Sequence dictionary for reference FASTA file"
        callingIntervalList: "Interval list of regions from which to call variants. Used for parallelization."
        knownVcfs: "Array of known indels VCF files"
        knownVcfIndexes: "Array of index files for known indels VCF files"
        dbSNP_vcf: "dbSNP VCF file"
        dbSNP_vcf_index: "Index file for dbSNP VCF file"
        scatterCount: "Number of intervals to scatter over"
        output_vcf_name: "Name of output VCF file"
    }

    input {
        File bam
        File bam_index
        File fasta
        File fasta_index
        File dict
        File callingIntervalList
        Array[File] knownVcfs
        Array[File] knownVcfIndexes
        File dbSNP_vcf
        File dbSNP_vcf_index
        Int scatterCount = 6
        String output_vcf_name = basename(bam, '.bam') + '.vcf.gz'
    }
    
    call picard.mark_duplicates {
        input: 
            bam = bam,
            create_bam = true
    }

    call gatk.split_n_cigar_reads {
        input:
            bam = select_first([mark_duplicates.duplicate_marked_bam, '']),
            bam_index = select_first([mark_duplicates.duplicate_marked_bam_index, '']),
            fasta = fasta,
            fasta_index = fasta_index,
            dict = dict,
            interval_list = callingIntervalList
    }

    call gatk.base_recalibrator {
        input:
            bam = split_n_cigar_reads.split_n_reads_bam,
            bam_index = split_n_cigar_reads.split_n_reads_bam_index,
            fasta = fasta,
            fasta_index = fasta_index,
            dict = dict,
            known_indels_sites_VCFs = knownVcfs,
            known_indels_sites_indices = knownVcfIndexes,
            dbSNP_vcf = dbSNP_vcf,
            dbSNP_vcf_index = dbSNP_vcf_index
    }

    call gatk.apply_bqsr {
        input:
            bam = split_n_cigar_reads.split_n_reads_bam,
            bam_index = split_n_cigar_reads.split_n_reads_bam_index,
            recalibration_report = base_recalibrator.recalibration_report
    }

    call picard.scatter_interval_list {
        input:
            interval_list = callingIntervalList,
            scatter_count = scatterCount
    }

    scatter (i in range(scatter_interval_list.interval_count)) {
        call gatk.haplotype_caller {
            input:
                bam = apply_bqsr.recalibrated_bam,
                bam_index = apply_bqsr.recalibrated_bam_index,
                fasta = fasta,
                fasta_index = fasta_index,
                dict = dict,
                interval_list = scatter_interval_list.out[i],
                dbSNP_vcf = dbSNP_vcf,
                dbSNP_vcf_index = dbSNP_vcf_index
        }
    }

    call picard.merge_vcfs {
        input:
            vcfs = haplotype_caller.vcf,
            vcfs_indexes = haplotype_caller.vcf_index,
            output_vcf_name = output_vcf_name
    }

    call gatk.variant_filtration {
        input: 
            vcf = merge_vcfs.output_vcf,
            vcf_index = merge_vcfs.output_vcf_index,
            fasta = fasta,
            fasta_index = fasta_index,
            dict = dict
    }

    output {
        File recalibrated_bam = apply_bqsr.recalibrated_bam
        File recalibrated_bam_index = apply_bqsr.recalibrated_bam_index
        File merged_vcf = merge_vcfs.output_vcf
        File merged_vcf_index = merge_vcfs.output_vcf_index
        File variant_filtered_vcf = variant_filtration.vcf_filtered
        File variant_filtered_vcf_index = variant_filtration.vcf_filtered_index
    }
}

