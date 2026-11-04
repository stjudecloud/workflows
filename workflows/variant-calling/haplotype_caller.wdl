version 1.3

import "../../tools/gatk4.wdl"
import "../../tools/picard.wdl"

workflow haplotype_caller {
    meta {
        description: "Workflow for calling germline variants using GATK's HaplotypeCaller."
        outputs: {
            raw_gvcf: "Raw gVCF file output by HaplotypeCaller",
            raw_gvcf_index: "Index file for the raw gVCF",
            raw_vcf: "Raw VCF file output by GenotypeGVCFs",
            raw_vcf_index: "Index file for the raw VCF",
            recalibrated_vcf: "Recalibrated VCF file after applying VQSR",
            recalibrated_vcf_index: "Index file for the recalibrated VCF",
            vcf_final: "Final VCF file after calculating genotype posteriors",
            vcf_final_index: "Index file for the final VCF",
        }
    }

    parameter_meta {
        bam: "Input BAM format file on which to perform germline variant calling"
        bam_index: "BAM index file corresponding to the input BAM"
        reference_fasta: "Reference genome in FASTA format"
        reference_fasta_index: "Index for FASTA format genome"
        reference_dict: "Dictionary file for FASTA format genome"
        dbSNP_vcf: "dbSNP VCF file"
        dbSNP_vcf_index: "dbSNP VCF index file"
        interval_list: "Interval list indicating regions in which to call variants"
        known_indels_sites_vcfs: "List of VCF files containing known indel sites to use for base quality score recalibration"
        known_indels_sites_indices: "List of VCF index files corresponding to the VCF files in `known_indels_sites_vcfs`"
        resources: {
            description: "List of resources to use for building the variant quality score recalibration model.",
            help: " Each resource should be a tuple containing the name of the resource, the path to the VCF file for the resource, and the path to the VCF index file for the resource.",
        }
        prefix: "Prefix for the output files."
    }

    input {
        File bam
        File bam_index
        File reference_fasta
        File reference_fasta_index
        File reference_dict
        #@ except: NamingConvention
        File dbSNP_vcf
        #@ except: NamingConvention
        File dbSNP_vcf_index
        File interval_list
        Array[File] known_indels_sites_vcfs
        Array[File] known_indels_sites_indices
        Array[Resource] resources
        String prefix = basename(bam, ".bam")
    }

    scatter (resource in resources) {
        call gatk4.resource_to_string {
            res = resource,
        }
    }

    call gatk4.base_recalibrator {
        bam,
        bam_index,
        fasta = reference_fasta,
        fasta_index = reference_fasta_index,
        dict = reference_dict,
        dbSNP_vcf,
        dbSNP_vcf_index,
        known_indels_sites_vcfs,
        known_indels_sites_indices,
        outfile_name = prefix + ".recalibration_report.txt",
    }

    call gatk4.apply_bqsr {
        bam,
        bam_index,
        recalibration_report = base_recalibrator.recalibration_report,
        prefix = prefix + ".recal",
    }

    call picard.scatter_interval_list {
        interval_list,
        scatter_count = 23,
    # break_bands_at_multiples_of = 10,
    }

    scatter (index in range(scatter_interval_list.interval_count)) {
        call gatk4.haplotype_caller {
            bam = apply_bqsr.recalibrated_bam,
            bam_index = apply_bqsr.recalibrated_bam_index,
            interval_list,
            fasta = reference_fasta,
            fasta_index = reference_fasta_index,
            dict = reference_dict,
            dbSNP_vcf,
            dbSNP_vcf_index,
            prefix,
            reference_confidence = RefConfidence.GVCF,
        }
    }

    call picard.merge_vcfs {
        vcfs = haplotype_caller.vcf,
        vcfs_indexes = haplotype_caller.vcf_index,
        output_vcf_name = prefix + ".vcf.gz",
    }

    call gatk4.genotype_gvcfs {
        gvcf = merge_vcfs.merged_vcf,
        gvcf_index = merge_vcfs.merged_vcf_index,
        fasta = reference_fasta,
        fasta_index = reference_fasta_index,
        dict = reference_dict,
        prefix,
    }

    call gatk4.variant_recalibrator {
        reference_fasta,
        reference_fasta_index,
        reference_dict,
        vcf = genotype_gvcfs.vcf,
        vcf_index = genotype_gvcfs.vcf_index,
        resources = resource_to_string.res_string,
        resource_vcfs = resource_to_string.res_vcf,
        resource_vcf_indices = resource_to_string.res_vcf_index,
    }

    call gatk4.apply_vqsr {
        reference_fasta,
        reference_fasta_index,
        reference_dict,
        vcf = genotype_gvcfs.vcf,
        vcf_index = genotype_gvcfs.vcf_index,
        recal_file = variant_recalibrator.recal_file,
        recal_file_index = variant_recalibrator.recal_index,
        tranches_file = variant_recalibrator.tranches_file,
        prefix = prefix + ".vqsr",
    }

    call gatk4.calculate_genotype_posteriors {
        vcf = apply_vqsr.vcf_recalibrated,
        vcf_index = apply_vqsr.vcf_recalibrated_index,
        supporting_vcf = dbSNP_vcf,
        supporting_vcf_index = dbSNP_vcf_index,
        prefix = prefix + ".posteriors",
    }

    output {
        File raw_gvcf = merge_vcfs.merged_vcf
        File raw_gvcf_index = merge_vcfs.merged_vcf_index
        File raw_vcf = genotype_gvcfs.vcf
        File raw_vcf_index = genotype_gvcfs.vcf_index
        File recalibrated_vcf = apply_vqsr.vcf_recalibrated
        File recalibrated_vcf_index = apply_vqsr.vcf_recalibrated_index
        File vcf_final = calculate_genotype_posteriors.vcf_posteriors
        File vcf_final_index = calculate_genotype_posteriors.vcf_posteriors_index
    }
}
