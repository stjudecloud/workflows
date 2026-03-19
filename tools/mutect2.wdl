version 1.3

task mutect2 {
    meta {
        description: "Run GATK Mutect2 somatic variant calling workflow"
        outputs: {
            mutect2_output: "Directory containing Mutect2 somatic variant calls",
            somatic_vcf: "VCF file with somatic variants called by Mutect2",
            somatic_vcf_index: "Index file for the Mutect2 somatic variants VCF",
            log_file: "Log file from the Mutect2 workflow execution",
        }
    }

    parameter_meta {
        reference_fasta: "Reference genome in FASTA format"
        reference_fasta_index: "Index file for the reference genome FASTA"
        normal_bam: "Input BAM file with aligned reads for normal sample"
        normal_bam_index: "Index file for the normal BAM file"
        tumor_bam: "Input BAM file with aligned reads for tumor sample"
        tumor_bam_index: "Index file for the tumor BAM file"
        germline_resource_vcf: "Optional VCF file with germline variants for Mutect2, recommended to be from gnomAD or similar population resource"
        germline_resource_vcf_index: "Index file for the germline resource VCF"
        panel_of_normals_vcf: "Optional VCF file with panel of normals for Mutect2, recommended to be generated from a large set of normal samples processed with Mutect2"
        panel_of_normals_vcf_index: "Index file for the panel of normals VCF"
        normal_sample_name: "Name of the normal sample"
        tumor_sample_name: "Name of the tumor sample"
        output_dir: "Directory to store Mutect2 output"
        threads: "Number of threads to use"
        modify_disk_size_gb: "Additional disk size in GB to allocate"
    }

    input {
        File reference_fasta
        File reference_fasta_index
        File normal_bam
        File normal_bam_index
        File tumor_bam
        File tumor_bam_index
        File? germline_resource_vcf
        File? germline_resource_vcf_index
        File? panel_of_normals_vcf
        File? panel_of_normals_vcf_index
        String normal_sample_name = basename(normal_bam, ".bam")
        String tumor_sample_name = basename(tumor_bam, ".bam")
        String output_dir = "mutect2_somatic_output"
        Int threads = 4
        Int modify_disk_size_gb = 0
    }

    Int disk_size_gb = ceil(size(reference_fasta, "GB") * 2)
        + ceil(size(normal_bam, "GB") * 2)
        + ceil(size(tumor_bam, "GB") * 2)
        + 20
        + modify_disk_size_gb

    String tumor = basename(tumor_bam)
    String normal = basename(normal_bam)

    command <<<
        set -euo pipefail

        ref_fasta=~{basename(reference_fasta, ".gz")}
        gunzip -c "~{reference_fasta}" > "$ref_fasta" \
            || ln -sf "~{reference_fasta}" "$ref_fasta"
        ln -sf "~{reference_fasta_index}" "$ref_fasta.fai"

        ln -s "~{tumor_bam}" "~{tumor}"
        ln -s "~{tumor_bam_index}" "~{tumor}.bai"
        ln -s "~{normal_bam}" "~{normal}"
        ln -s "~{normal_bam_index}" "~{normal}.bai"

        gatk Mutect2 \
            -R "$ref_fasta" \
            -I "~{tumor}" \
            -I "~{normal}" \
            -normal "~{normal_sample_name}" \
            -tumor "~{tumor_sample_name}" \
            ~{if defined(germline_resource_vcf) then "-germline-resource '" + germline_resource_vcf + "'" else ""} \
            ~{if defined(panel_of_normals_vcf) then "-panel-of-normals '" + panel_of_normals_vcf + "'" else ""} \
            -O "~{output_dir}/somatic_variants.vcf.gz" \
            --native-pair-hmm-threads "~{threads}"

        rm -rf "$ref_fasta" "$ref_fasta.fai" "~{tumor}" "~{tumor}.bai" "~{normal}" "~{normal}.bai"
     >>>

     output {
         Directory mutect2_output = output_dir
         File somatic_vcf = "~{output_dir}/somatic_variants.vcf.gz"
         File somatic_vcf_index = "~{output_dir}/somatic_variants.vcf.gz.tbi"
         File log_file = "~{output_dir}/gatk_mutect2.log"
    }

    requirements {
        container: "quay.io/biocontainers/gatk4:4.6.2.0--py310hdfd78af_1"
        cpu: threads
        memory: "25 GB"
        disks: "~{disk_size_gb} GB"
    }
}
