version 1.2

task clair3 {
    meta {
        description: "Run Clair3 variant caller for small variants using deep neural networks"
        outputs: {
            pileup_vcf: "VCF file with variants called using pileup model",
            full_alignment_vcf: "VCF file with variants called using full-alignment model",
            merged_vcf: "Final merged VCF file with variants from both models",
        }
    }

    parameter_meta {
        reference_fasta: "Reference genome in FASTA format"
        bam: "Input BAM file with aligned reads"
        model: "Pre-trained Clair3 model to use for variant calling"
        bed_regions: "Optional BED file specifying regions to call variants in"
        vcf_candidates: "Optional VCF file with candidate variants to consider"
        output_dir: "Directory to store Clair3 output"
        platform: {
            description: "Sequencing platform used to generate the reads",
            choices: [
                "ont",
                "hifi",
                "ilmn",
            ],
        }
        all_contigs: "Boolean indicating whether to include all contigs in variant calling. If false only chr{1..22,X,Y} are called."
        print_ref_calls: "Boolean indicating whether to print reference calls in the output VCF"
        gvcf: "Boolean indicating whether to output gVCF format"
        threads: "Number of threads to use"
        modify_disk_size_gb: "Additional disk size in GB to allocate"
    }

    input {
        File reference_fasta
        File bam
        File model
        File? bed_regions
        File? vcf_candidates
        String output_dir = "clair3_output"
        String platform = "ilmn"
        Boolean all_contigs = false
        Boolean print_ref_calls = false
        Boolean gvcf = false
        Int threads = 4
        Int modify_disk_size_gb = 0
    }

    Int disk_size_gb = ceil(size(reference_fasta, "GiB") * 2)
        + ceil(size(bam, "GiB"))
        + 20
        + modify_disk_size_gb

    command <<<
        set -euo pipefail

        ref_fasta=~{basename(reference_fasta, ".gz")}
        gunzip -c "~{reference_fasta}" > "$ref_fasta" \
            || ln -sf "~{reference_fasta}" "$ref_fasta"

        ./run_clair3.sh \
            --bam_fn="~{bam}" \
            --ref_fn="$ref_fasta" \
            --threads="~{threads}" \
            --platform="~{platform}" \
            --model_path="~{model}" \
            --output="~{output_dir}" \
            ~{if all_contigs then "--include_all_ctgs" else ""} \
            ~{if print_ref_calls then "--print_ref_calls" else ""} \
            ~{if defined(bed_regions) then "--bed_fn='~{bed_regions}'" else ""} \
            ~{if defined(vcf_candidates) then "--vcf_fn='~{vcf_candidates}'" else ""} \
            ~{if gvcf then "--gvcf" else ""}

        rm -rf "$ref_fasta"
    >>>

    output {
        File pileup_vcf = "~{output_dir}/pileup.vcf.gz"
        File full_alignment_vcf = "~{output_dir}/full_alignment.vcf.gz"
        File merged_vcf = "~{output_dir}/merge_output.vcf.gz"
    }

    requirements {
        container: "quay.io/biocontainers/clair3:1.2.0--py310h779eee5_0"
        cpu: threads
        memory: "16 GB"
        disks: "~{disk_size_gb} GB"
    }
}
