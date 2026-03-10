version 1.2

task germline {
    meta {
        description: "Run Octopus individual germline variant caller"
        outputs: {
            output_vcf: "VCF file with called germline variants"
        }
    }

    parameter_meta {
        reference_fasta: "Reference genome in FASTA format"
        reference_fasta_index: "Index file for the reference genome FASTA"
        bam: "Input BAM file with aligned reads"
        forest: "Pre-trained random forest model for variant filtering"
        output_vcf: "Output VCF file with called variants"
        threads: "Number of threads to use"
        modify_disk_size_gb: "Additional disk size in GB to allocate"
    }

    input {
        File reference_fasta
        File reference_fasta_index
        File bam
        File forest
        String output_vcf = "octopus_germline_output.vcf"
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
        ln -sf "~{reference_fasta_index}" "$ref_fasta.fai"

        octopus \
            -R "$ref_fasta" \
            -I "~{bam}" \
            -o "~{output_vcf}" \
            --forest "~{forest}" \
            --threads "~{threads}"

        rm -rf "$ref_fasta" "$ref_fasta.fai"
    >>>

    output {
        File output_vcf = "~{output_vcf}"
    }

    requirements {
        container: "quay.io/biocontainers/octopus:0.7.4--ha3c1580_2"
        cpu: threads
        memory: "30 GB"
        disks: "~{disk_size_gb} GB"
}

task somatic {
    meta {
        description: "Run Octopus individual somatic variant caller"
        outputs: {
            output_vcf: "VCF file with called somatic variants"
        }
    }

    parameter_meta {
        reference_fasta: "Reference genome in FASTA format"
        reference_fasta_index: "Index file for the reference genome FASTA"
        normal_bam: "Input BAM file with aligned reads from normal sample"
        tumor_bam: "Input BAM file with aligned reads from tumor sample"
        forest: "Pre-trained random forest model for variant filtering"
        output_vcf: "Output VCF file with called variants"
        threads: "Number of threads to use"
        modify_disk_size_gb: "Additional disk size in GB to allocate"
    }

    input {
        File reference_fasta
        File reference_fasta_index
        File normal_bam
        File tumor_bam
        File forest
        String output_vcf = "octopus_germline_output.vcf"
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
        ln -sf "~{reference_fasta_index}" "$ref_fasta.fai"

        # TODO: I have no idea what `--normal-sample` is supposed to contain. The documentation is unhelpful.
        octopus \
            -R "$ref_fasta" \
            -I "~{normal_bam}" "~{tumor_bam}" \
            --normal-sample "~{basename(normal_bam, ".bam")}" \
            -o "~{output_vcf}" \
            --forest "~{forest}" \
            --threads "~{threads}"

        rm -rf "$ref_fasta"
    >>>

    output {
        File output_vcf = "~{output_vcf}"
    }

    requirements {
        container: "quay.io/biocontainers/octopus:0.7.4--ha3c1580_2"
        cpu: threads
        memory: "30 GB"
        disks: "~{disk_size_gb} GB"
}