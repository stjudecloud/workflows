version 1.2

task deepsomatic {
    meta {
        description: "Call variants using DeepSomatic"
        outputs: {
            vcf_output: "VCF file containing called somatic variants",
            gvcf_output: "gVCF file containing called somatic variants",
            runtime: "Optional HTML report of runtime metrics",
            vcf_stats: "Optional HTML report of VCF statistics",
        }
    }

    parameter_meta {
        reference_fasta: "Reference genome in FASTA format"
        reference_fasta_index: "Index file for the reference genome FASTA"
        tumor_bam: "Input BAM file with aligned reads for tumor sample"
        normal_bam: "Input BAM file with aligned reads for normal sample"
        output_prefix: "Prefix for output VCF and gVCF files"
        tumor_sample_name: "Sample name for the tumor sample"
        normal_sample_name: "Sample name for the normal sample"
        model_type: {
            description: "Type of model to use for variant calling",
            choices: [
                "WGS",
                "WES",
                "PACBIO",
                "ONT",
                "FFPE_WGS",
                "FFPE_WES",
                "FFPE_WGS_TUMOR_ONLY",
                "FFPE_WES_TUMOR_ONLY",
                "WGS_TUMOR_ONLY",
                "WES_TUMOR_ONLY",
                "PACBIO_TUMOR_ONLY",
                "ONT_TUMOR_ONLY",
            ],
        }
        runtime_report: "Output make_examples_somatic runtime metrics and create a visual runtime report using runtime_by_region_vis."
        vcf_stats_report: "Output a visual report (HTML) of statistics about the output VCF."
        threads: "Number of threads to use"
        modify_disk_size_gb: "Additional disk size in GB to allocate"
    }

    input {
        File reference_fasta
        File reference_fasta_index
        File tumor_bam
        File normal_bam
        String output_prefix = "deepsomatic_output"
        String tumor_sample_name = "tumor"
        String normal_sample_name = "normal"
        String model_type = "WGS"
        Boolean runtime_report = false
        Boolean vcf_stats_report = false
        Int threads = 8
        Int modify_disk_size_gb = 0
    }

    Int disk_size_gb = ceil(size(reference_fasta, "GiB") * 2)
        + ceil(size(tumor_bam, "GiB"))
        + ceil(size(normal_bam, "GiB"))
        + 50
        + modify_disk_size_gb

    command <<<
        set -euo pipefail

        ref_fasta=~{basename(reference_fasta, ".gz")}
        gunzip -c "~{reference_fasta}" > "$ref_fasta" \
            || ln -sf "~{reference_fasta}" "$ref_fasta"
        ln -sf "~{reference_fasta_index}" "$ref_fasta.fai"

        run_deepsomatic \
            --model_type="~{model_type}" \
            --ref="$ref_fasta" \
            --tumor_bam="~{tumor_bam}" \
            --normal_bam="~{normal_bam}" \
            --output_vcf="~{output_prefix}.vcf.gz" \
            --output_gvcf="~{output_prefix}.g.vcf.gz" \
            --tumor_sample_name="~{tumor_sample_name}" \
            --normal_sample_name="~{normal_sample_name}" \
            --num_shards="~{threads}" \
            --logging_dir="logs" \
            --intermediate_results_dir="intermediate_results" \
            ~{if runtime_report then "--runtime_report" else ""} \
            ~{if vcf_stats_report then "--vcf_stats_report" else ""}


        rm -rf "$ref_fasta"

    >>>

    output {
        File vcf_output = "~{output_prefix}.vcf.gz"
        File gvcf_output = "~{output_prefix}.g.vcf.gz"
        File? runtime = "logs/runtime_by_region_vis.html"
        File? vcf_stats = "logs/vcf_stats_report.html"
    }

    requirements {
        container: "google/deepsomatic:1.9.0-gpu"
        cpu: threads
        memory: "32 GB"
        disks: "~{disk_size_gb} GB"
        gpu: true
    }

    hints {
        gpu: 1
    }
}

task deepvariant {
    meta {
        description: "Call variants using DeepVariant"
        outputs: {
            vcf_output: "VCF file containing called variants",
            gvcf_output: "gVCF file containing called variants",
            runtime: "Optional HTML report of runtime metrics",
            vcf_stats: "Optional HTML report of VCF statistics",
        }
    }

    parameter_meta {
        reference_fasta: "Reference genome in FASTA format"
        reference_fasta_index: "Index file for the reference genome FASTA"
        bam: "Input BAM file with aligned reads for sample"
        bam_index: "Index for input BAM file"
        haploid_chromosomes: "List of chromosomes to be treated as haploid during variant calling"
        output_prefix: "Prefix for output VCF and gVCF files"
        model_type: {
            description: "Type of model to use for variant calling",
            choices: [
                "WGS",
                "WES",
                "PACBIO",
                "ONT",
                "FFPE_WGS",
                "FFPE_WES",
                "FFPE_WGS_TUMOR_ONLY",
                "FFPE_WES_TUMOR_ONLY",
                "WGS_TUMOR_ONLY",
                "WES_TUMOR_ONLY",
                "PACBIO_TUMOR_ONLY",
                "ONT_TUMOR_ONLY",
            ],
        }
        runtime_report: "Output make_examples_somatic runtime metrics and create a visual runtime report using runtime_by_region_vis."
        vcf_stats_report: "Output a visual report (HTML) of statistics about the output VCF."
        threads: "Number of threads to use"
        modify_disk_size_gb: "Additional disk size in GB to allocate"
    }

    input {
        File reference_fasta
        File reference_fasta_index
        File bam
        File bam_index
        Array[String] haploid_chromosomes = ["chrX", "chrY"]
        String output_prefix = "deepsomatic_output"
        String model_type = "WGS"
        Boolean runtime_report = false
        Boolean vcf_stats_report = false
        Int threads = 8
        Int modify_disk_size_gb = 0
    }

    Int disk_size_gb = ceil(size(reference_fasta, "GiB") * 2)
        + ceil(size(bam, "GiB"))
        + 50
        + modify_disk_size_gb

    command <<<
        set -euo pipefail

        ref_fasta=~{basename(reference_fasta, ".gz")}
        gunzip -c "~{reference_fasta}" > "$ref_fasta" \
            || ln -sf "~{reference_fasta}" "$ref_fasta"
        ln -sf "~{reference_fasta_index}" "$ref_fasta.fai"

        run_deepvariant \
            --model_type="~{model_type}" \
            --ref="$ref_fasta" \
            --reads="~{bam}" \
            --output_vcf="~{output_prefix}.vcf.gz" \
            --output_gvcf="~{output_prefix}.g.vcf.gz" \
            --num_shards="~{threads}" \
            --logging_dir="logs" \
            --intermediate_results_dir="intermediate_results" \
            ~{if runtime_report then "--runtime_report" else ""} \
            ~{if vcf_stats_report then "--vcf_stats_report" else ""} \
            --haploid_contigs="~{sep(",", haploid_chromosomes)}"

        rm -rf "$ref_fasta"
    >>>

    output {
        File vcf_output = "~{output_prefix}.vcf.gz"
        File gvcf_output = "~{output_prefix}.g.vcf.gz"
        File? runtime = "logs/runtime_by_region_vis.html"
        File? vcf_stats = "logs/vcf_stats_report.html"
    }

    requirements {
        container: "google/deepvariant:1.9.0-gpu"
        cpu: threads
        memory: "32 GB"
        disks: "~{disk_size_gb} GB"
        gpu: true
    }

    hints {
        gpu: 1
    }
}
