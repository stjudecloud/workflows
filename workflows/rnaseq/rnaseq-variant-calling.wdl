version 1.1

import "../../tools/gatk4.wdl" as gatk
import "../../tools/picard.wdl"

workflow rnaseq_variant_calling {
    meta {
        name: "RNA-Seq Variant Calling"
        description: "Call short germline variants from RNA-Seq data. Produces a VCF file of variants. Based on GATK RNA-Seq short variant calling best practices pipeline."
        category: "Variant Calling"
        outputs: {
            recalibrated_bam: "BAM that has undergone recalibration of base quality scores",
            recalibrated_bam_index: "Index file for recalibrated BAM file",
            variant_filtered_vcf: "VCF file after variant filters have been applied",
            variant_filtered_vcf_index: "Index for filtered variant VCF file",
        }
        allowNestedInputs: true
    }

    parameter_meta {
        bam: "BAM file of aligned RNA-Seq reads"
        bam_index: "Index file for BAM file"
        fasta: "Reference FASTA file"
        fasta_index: "Index file for reference FASTA file"
        dict: "Sequence dictionary for reference FASTA file"
        calling_interval_list: "Interval list of regions from which to call variants. Used for parallelization."
        dbSNP_vcf: "dbSNP VCF file"
        dbSNP_vcf_index: "Index file for dbSNP VCF file"
        known_vcfs: "Array of known indels VCF files"
        known_vcf_indexes: "Array of index files for known indels VCF files"
        prefix: "Prefix for the output files."
        bam_is_dup_marked: "Whether the input BAM file has duplicates marked."
        scatter_count: "Number of intervals to scatter over. This should typically be set to 5-20. Higher values will increase parallelism and speed up the workflow, but increase overhead in provisioning resources."
    }

    input {
        File bam
        File bam_index
        File fasta
        File fasta_index
        File dict
        File calling_interval_list
        #@ except: SnakeCase
        File dbSNP_vcf
        #@ except: SnakeCase
        File dbSNP_vcf_index
        Array[File] known_vcfs
        Array[File] known_vcf_indexes
        String prefix = basename(bam, ".bam")
        Boolean bam_is_dup_marked = false
        Int scatter_count = 6
    }

    if (!bam_is_dup_marked){
        call picard.mark_duplicates { input:
            bam,
            create_bam = true,
        }
    }

    call gatk.split_n_cigar_reads { input:
        bam = select_first([mark_duplicates.duplicate_marked_bam, bam]),
        bam_index = select_first([mark_duplicates.duplicate_marked_bam_index, bam_index]),
        fasta,
        fasta_index,
        dict,
    }

    call gatk.base_recalibrator { input:
        bam = split_n_cigar_reads.split_n_reads_bam,
        bam_index = split_n_cigar_reads.split_n_reads_bam_index,
        fasta,
        fasta_index,
        dict,
        known_indels_sites_vcfs = known_vcfs,
        known_indels_sites_indices = known_vcf_indexes,
        dbSNP_vcf,
        dbSNP_vcf_index,
    }

    call gatk.apply_bqsr { input:
        bam = split_n_cigar_reads.split_n_reads_bam,
        bam_index = split_n_cigar_reads.split_n_reads_bam_index,
        recalibration_report = base_recalibrator.recalibration_report,
        prefix,
    }

    call picard.scatter_interval_list { input:
        interval_list = calling_interval_list,
        scatter_count,
    }

    scatter (list in scatter_interval_list.interval_lists_scatter) {
        call gatk.haplotype_caller { input:
            bam = apply_bqsr.recalibrated_bam,
            bam_index = apply_bqsr.recalibrated_bam_index,
            fasta,
            fasta_index,
            dict,
            interval_list = list,
            dbSNP_vcf,
            dbSNP_vcf_index,
        }
    }

    call picard.merge_vcfs { input:
        vcfs = haplotype_caller.vcf,
        vcfs_indexes = haplotype_caller.vcf_index,
        output_vcf_name = "~{prefix}.vcf.gz",
    }

    call gatk.variant_filtration { input:
        vcf = merge_vcfs.merged_vcf,
        vcf_index = merge_vcfs.merged_vcf_index,
        fasta,
        fasta_index,
        dict,
        prefix,
    }

    output {
        File recalibrated_bam = apply_bqsr.recalibrated_bam
        File recalibrated_bam_index = apply_bqsr.recalibrated_bam_index
        File variant_filtered_vcf = variant_filtration.vcf_filtered
        File variant_filtered_vcf_index = variant_filtration.vcf_filtered_index
    }
}
