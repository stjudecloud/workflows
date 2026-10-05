version 1.3

import "../../tools/mutect2.wdl" as mutect2_task

workflow mutect2 {
    meta {
        description: "Workflow for calling somatic variants using GATK Mutect2 and filtering with FilterMutectCalls"
        outputs: {
            filtered_somatic_vcf: "VCF file with filtered somatic variants from Mutect2",
            filtered_somatic_vcf_index: "Index file for the filtered Mutect2 somatic variants VCF",
        }
    }

    parameter_meta {
        reference_fasta: "Reference genome in FASTA format"
        reference_fasta_index: "Index file for the reference genome FASTA"
        reference_fasta_dict: "Dictionary file for the reference genome FASTA"
        normal_bam: "Input BAM file with aligned reads for normal sample"
        normal_bam_index: "Index file for the normal BAM file"
        tumor_bam: "Input BAM file with aligned reads for tumor sample"
        tumor_bam_index: "Index file for the tumor BAM file"
        variant_vcf: "VCF file with variants and allele frequencies to summarize pileups over."
        variant_vcf_index: "Index file for the variant VCF"
        intervals: "One or more genomic intervals over which to operate. Often the same file as `variants`"
        intervals_index: "Index file for the intervals file"
        germline_resource_vcf: "Optional VCF file with germline variants for Mutect2, recommended to be from gnomAD or similar population resource"
        germline_resource_vcf_index: "Index file for the germline resource VCF"
        panel_of_normals_vcf: "Optional VCF file with panel of normals for Mutect2, recommended to be generated from a large set of normal samples processed with Mutect2"
        panel_of_normals_vcf_index: "Index file for the panel of normals VCF"
        normal_sample_name: "Name of the normal sample"
        tumor_sample_name: "Name of the tumor sample"
        output_prefix: "Prefix for output files. The extensions '.vcf.gz' and '.vcf.gz.tbi' will be added."
    }

    input {
        File reference_fasta
        File reference_fasta_index
        File reference_fasta_dict
        File normal_bam
        File normal_bam_index
        File tumor_bam
        File tumor_bam_index
        File variant_vcf
        File variant_vcf_index
        File intervals
        File intervals_index
        File? germline_resource_vcf
        File? germline_resource_vcf_index
        File? panel_of_normals_vcf
        File? panel_of_normals_vcf_index
        String normal_sample_name = basename(normal_bam, ".bam")
        String tumor_sample_name = basename(tumor_bam, ".bam")
        String output_prefix = basename(tumor_bam, ".bam") + "_v_" + basename(normal_bam, ".bam"
        )
    }

    call mutect2_task.mutect2 {
        reference_fasta,
        reference_fasta_index,
        reference_fasta_dict,
        normal_bam,
        normal_bam_index,
        tumor_bam,
        tumor_bam_index,
        normal_sample_name,
        tumor_sample_name,
        output_prefix = output_prefix + "_unfiltered",
        germline_resource_vcf,
        germline_resource_vcf_index,
        panel_of_normals_vcf,
        panel_of_normals_vcf_index,
    }

    call mutect2_task.get_pileup_summaries as get_tumor_pileups {
        bam = tumor_bam,
        bam_index = tumor_bam_index,
        intervals,
        intervals_index,
        variants = variant_vcf,
        variants_index = variant_vcf_index,
        prefix = output_prefix + "_tumor",
        reference_fasta,
        reference_fasta_index,
        reference_fasta_dict,
    }

    call mutect2_task.get_pileup_summaries as get_normal_pileups {
        bam = normal_bam,
        bam_index = normal_bam_index,
        intervals,
        intervals_index,
        variants = variant_vcf,
        variants_index = variant_vcf_index,
        prefix = output_prefix + "_normal",
        reference_fasta,
        reference_fasta_index,
        reference_fasta_dict,
    }

    call mutect2_task.calculate_contamination {
        tumor_pileups = get_tumor_pileups.pileup_summaries,
        normal_pileups = get_normal_pileups.pileup_summaries,
        prefix = output_prefix,
    }

    call mutect2_task.filter_mutect {
        unfiltered_somatic_vcf = mutect2.somatic_vcf,
        unfiltered_somatic_vcf_index = mutect2.somatic_vcf_index,
        unfiltered_somatic_vcf_stats = mutect2.stats,
        reference_fasta,
        reference_fasta_index,
        reference_fasta_dict,
        contamination_table = calculate_contamination.contamination_table,
        maf_segments = calculate_contamination.maf_segments,
        prefix = output_prefix,
    }

    output {
        File filtered_somatic_vcf = filter_mutect.filtered_somatic_vcf
        File filtered_somatic_vcf_index = filter_mutect.filtered_somatic_vcf_index
    }
}

