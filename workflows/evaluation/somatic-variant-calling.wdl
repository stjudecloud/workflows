version 1.3

import "../../data_structures/read_group.wdl"
import "../../tools/clair.wdl"
import "../../tools/deepvariant.wdl"
import "../../tools/manta.wdl"
import "../../tools/mutect2.wdl"
import "../../tools/strelka.wdl"

workflow somatic_variant_calling {
    meta {
        description: "Workflow for calling somatic variants using multiple tools"
        outputs: {
            clairs_vcf: "VCF file with called somatic variants from Clair",
            deepsomatic_vcf: "VCF file with called somatic variants from DeepSomatic",
            manta_output: "Directory containing Manta variant calls",
            strelka_output: "Directory containing Strelka somatic variant calls",
            mutect2_vcf: "VCF file with filtered somatic variants from Mutect2",
            mutect2_vcf_index: "Index file for the Mutect2 filtered somatic variants VCF",
        }
        allowNestedInputs: true
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
        deepsomatic_model_type: "Type of model to use for DeepSomatic variant calling"
        exome: "Whether the data is from exome sequencing, which will adjust filtering in some tools accordingly"
        run_deepsomatic: "Whether to run DeepSomatic for variant calling"
        run_strelka: "Whether to run Strelka for variant calling"
        run_clair_s: "Whether to run Clair-S for variant calling"
        run_mutect2: "Whether to run Mutect2 workflow for variant calling"
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
        ModelType deepsomatic_model_type = ModelType.WGS
        Boolean exome = false
        Boolean run_deepsomatic = true
        Boolean run_strelka = true
        Boolean run_clair_s = true
        Boolean run_mutect2 = true
    }

    call read_group.get_read_groups as normal_read_groups {
        bam = normal_bam,
    }

    call read_group.get_read_groups as tumor_read_groups {
        bam = tumor_bam,
    }

    if (run_deepsomatic) {
        call deepvariant.deepsomatic {
            reference_fasta,
            reference_fasta_index,
            tumor_bam,
            tumor_bam_index,
            normal_bam,
            normal_bam_index,
            output_prefix = "deepsomatic_output",
            tumor_sample_name = "tumor",
            normal_sample_name = "normal",
            model_type = deepsomatic_model_type,
            runtime_report = true,
            vcf_stats_report = true,
        }
    }

    if (run_strelka) {
        call manta.manta_somatic {
            reference_fasta,
            reference_fasta_index,
            normal_bam,
            normal_bam_index,
            tumor_bam,
            tumor_bam_index,
            output_dir = "manta_output",
            exome,
        }

        call strelka.somatic {
            reference_fasta,
            reference_fasta_index,
            normal_bam,
            normal_bam_index,
            tumor_bam,
            tumor_bam_index,
            indel_candidates = manta_somatic.indel_candidates,
            indel_candidates_index = manta_somatic.indel_candidates_index,
            output_dir = "strelka_output",
            exome,
        }
    }

    if (run_clair_s) {
        call clair.clair_s {
            reference_fasta,
            reference_fasta_index,
            normal_bam,
            normal_bam_index,
            tumor_bam,
            tumor_bam_index,
            platform = ClairSPlatform.ilmn,
            prefix = "~{basename(tumor_bam, ".bam")}_vs_~{basename(normal_bam, ".bam")}",
            sample_name = "~{basename(tumor_bam, ".bam")}",
        }
    }

    if (run_mutect2) {
        call mutect2.mutect2_wf {
            reference_fasta,
            reference_fasta_index,
            reference_fasta_dict,
            normal_bam,
            normal_bam_index,
            tumor_bam,
            tumor_bam_index,
            normal_sample_name = select_first([normal_read_groups.read_groups[0].SM, "normal"]
            ),
            tumor_sample_name = select_first([tumor_read_groups.read_groups[0].SM, "tumor"]),
            variant_vcf,
            variant_vcf_index,
            intervals,
            intervals_index,
        }
    }

    output {
        File? clairs_vcf = clair_s.vcf
        File? deepsomatic_vcf = deepsomatic.vcf_output
        Directory? manta_output = manta_somatic.manta_output
        Directory? strelka_output = somatic.strelka_output
        File? mutect2_vcf = mutect2_wf.filtered_somatic_vcf
        File? mutect2_vcf_index = mutect2_wf.filtered_somatic_vcf_index
    }
}
