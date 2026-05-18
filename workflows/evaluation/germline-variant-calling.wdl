version 1.3

import "../../tools/clair.wdl"
import "../../tools/deepvariant.wdl"
import "../../tools/gatk4.wdl"
import "../../tools/manta.wdl"
import "../../tools/ngsep.wdl"
import "../../tools/strelka.wdl"

workflow variant_calling {
    meta {
        description: "Runs a series of variant calling tools on a processed BAM file."
        help: "The workflow is designed to be flexible and allow users to choose which variant callers to run based on their specific needs and preferences."
        outputs: {
            clair_vcf: "VCF file produced by Clair",
            clair_full_vcf: "VCF file produced by Clair containing all variants before filtering",
            clair_merged_vcf: "VCF file produced by Clair after merging pileup and full-alignment calls",
            deepvariant_vcf: "VCF file produced by DeepVariant",
            deepvariant_gvcf: "gVCF file produced by DeepVariant",
            deepvariant_runtime_report: "Optional HTML report of runtime metrics from DeepVariant",
            deepvariant_vcf_stats: "Optional HTML report of VCF statistics from DeepVariant",
            manta_output: "Directory containing Manta structural variant calls and associated files",
            manta_log: "Log file from the Manta workflow execution",
            ngsep_vcf: "VCF file produced by NGSEP",
            strelka_output: "Directory containing Strelka somatic variant calls and associated files",
            strelka_log: "Log file from the Strelka workflow execution",
            haplotype_caller_vcf: "VCF file output by GATK HaplotypeCaller after VQSR and genotype posterior calculation",
            haplotype_caller_vcf_index: "Index file for the HaplotypeCaller VCF output",
        }
        allowNestedInputs: true
    }

    parameter_meta {
        processed_bam: "Input BAM format file that has been processed and is ready for variant calling"
        processed_bam_index: "BAI index file associated with the input BAM file"
        reference_genome: "Reference genome in FASTA format that corresponds to the alignments in the input BAM file"
        reference_genome_index: "Index file for the reference genome FASTA, typically with a .fai extension"
        reference_genome_dictionary: "Dictionary file for FASTA format genome"
        dbSNP_vcf: "dbSNP VCF file"
        dbSNP_vcf_index: "dbSNP VCF index file"
        interval_list: {
            description: "Interval list indicating regions in which to call variants",
            external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360035531852-Intervals-and-interval-lists",
        }
        clair3_model: "Pre-trained Clair3 model file to use for variant calling with Clair"
        known_indels_sites_vcfs: "Optional array of VCF files containing known indel sites to be used in variant calling and filtering"
        known_indels_sites_indices: "Optional array of index files corresponding to the known indel site VCF files"
        resources: "Optional array of additional resource files to be used in variant calling and filtering, such as additional VCFs or BED files"
        run_clair: "Whether to run Clair for variant calling"
        run_deepvariant: "Whether to run DeepVariant for variant calling"
        run_haplotype_caller: "Whether to run GATK's Haplotype Caller for variant calling"
        run_manta: "Whether to run Manta for structural variant calling"
        run_ngsep: "Whether to run NGSEP for variant calling"
        run_strelka: "Whether to run Strelka for variant calling"
    }

    input {
        File processed_bam
        File processed_bam_index
        File reference_genome
        File reference_genome_index
        File reference_genome_dictionary
        #@ except: SnakeCase
        File dbSNP_vcf
        #@ except: SnakeCase
        File dbSNP_vcf_index
        File interval_list
        Directory clair3_model
        Array[File] known_indels_sites_vcfs = []
        Array[File] known_indels_sites_indices = []
        Array[Resource] resources = []
        Boolean run_clair = true
        Boolean run_deepvariant = true
        Boolean run_haplotype_caller = true
        Boolean run_manta = true
        Boolean run_ngsep = true
        Boolean run_strelka = true
    }

    if (run_clair) {
        call clair.clair3 {
            bam = processed_bam,
            bam_index = processed_bam_index,
            reference_fasta = reference_genome,
            reference_fasta_index = reference_genome_index,
            model = clair3_model,
        }
    }

    if (run_deepvariant) {
        call deepvariant.deepvariant {
            bam = processed_bam,
            bam_index = processed_bam_index,
            reference_fasta = reference_genome,
            reference_fasta_index = reference_genome_index,
        }
    }

    if (run_manta) {
        call manta.manta_germline {
            bam = processed_bam,
            bam_index = processed_bam_index,
            reference_fasta = reference_genome,
            reference_fasta_index = reference_genome_index,
        }
    }

    if (run_ngsep) {
        call ngsep.germline_variant as ngsep_germline_variant {
            bam = processed_bam,
            reference_fasta = reference_genome,
        }
    }

    if (run_strelka) {
        call strelka.germline as strelka_germline {
            bam = processed_bam,
            bam_index = processed_bam_index,
            reference_fasta = reference_genome,
            reference_fasta_index = reference_genome_index,
        }
    }

    if (run_haplotype_caller) {
        call gatk4.germline_variant_calling_wf {
            bam = processed_bam,
            bam_index = processed_bam_index,
            interval_list,
            reference_fasta = reference_genome,
            reference_fasta_index = reference_genome_index,
            reference_dict = reference_genome_dictionary,
            dbSNP_vcf,
            dbSNP_vcf_index,
            known_indels_sites_vcfs,
            known_indels_sites_indices,
            resources,
        }
    }

    output {
        File? clair_vcf = clair3.pileup_vcf
        File? clair_full_vcf = clair3.full_alignment_vcf
        File? clair_merged_vcf = clair3.merged_vcf
        File? deepvariant_vcf = deepvariant.vcf_output
        File? deepvariant_gvcf = deepvariant.gvcf_output
        File? deepvariant_runtime_report = deepvariant.runtime_report
        File? deepvariant_vcf_stats = deepvariant.vcf_stats
        Directory? manta_output = manta_germline.manta_output
        File? manta_log = manta_germline.log_file
        Array[File]? ngsep_vcf = ngsep_germline_variant.vcf_output
        Directory? strelka_output = strelka_germline.strelka_output
        File? strelka_log = strelka_germline.log_file
        File? haplotype_caller_vcf = germline_variant_calling_wf.vcf_final
        File? haplotype_caller_vcf_index = germline_variant_calling_wf.vcf_final_index
    }
}
