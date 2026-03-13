version 1.2

task manta_germline {
    meta {
        description: "Run Manta structural variant and indel caller"
        outputs: {
            manta_output: "Directory containing Manta variant calls",
            log_file: "Log file from the Manta workflow execution",
        }
    }

    parameter_meta {
        reference_fasta: "Reference genome in FASTA format"
        reference_fasta_index: "Index file for the reference genome FASTA"
        bam: "Input BAM file with aligned reads"
        bam_index: "Index file for the input BAM file"
        calling_regions_bed: "Optional BED file specifying regions to call variants"
        output_dir: "Directory to store Manta output"
        threads: "Number of threads to use"
        modify_disk_size_gb: "Additional disk size in GB to allocate"
    }

    input {
        File reference_fasta
        File reference_fasta_index
        File bam
        File bam_index
        File? calling_regions_bed
        String output_dir = "manta_output"
        Int threads = 4
        Int modify_disk_size_gb = 0
    }

    Int disk_size_gb = ceil(size(reference_fasta, "GB") * 2)
        + ceil(size(bam, "GB"))
        + 20
        + modify_disk_size_gb

    String filename = basename(bam)

    command <<<
        set -euo pipefail

        ref_fasta=~{basename(reference_fasta, ".gz")}
        gunzip -c "~{reference_fasta}" > "$ref_fasta" \
            || ln -sf "~{reference_fasta}" "$ref_fasta"
        ln -sf "~{reference_fasta_index}" "$ref_fasta.fai"

        ln -s "~{bam}" "~{filename}"
        ln -s "~{bam_index}" "~{filename}.bai"

        configManta.py \
            --bam "~{filename}" \
            --referenceFasta "$ref_fasta" \
            ~{if defined(calling_regions_bed) then "--callRegions '" + calling_regions_bed + "'" else ""} \
            --runDir "~{output_dir}"

        "~{output_dir}/runWorkflow.py" -j "~{threads}"

        rm -rf "$ref_fasta" "$ref_fasta.fai" "~{filename}" "~{filename}.bai"
    >>>

    output {
        Directory manta_output = output_dir
        File log_file = "~{output_dir}/workspace/pyflow.data/logs/pyflow_log.txt"
    }

    requirements {
        container: "quay.io/biocontainers/manta:1.6.0--py27h9948957_6"
        cpu: threads
        memory: "25 GB"
        disks: "~{disk_size_gb} GB"
    }
}

task manta_somatic {
    meta {
        description: "Run Manta structural variant and indel caller in somatic mode"
        outputs: {
            manta_output: "Directory containing Manta variant calls",
            indel_candidates: "VCF file with candidate small indels",
            indel_candidates_index: "Index file for the candidate small indels VCF",
            log_file: "Log file from the Manta workflow execution",
        }
    }

    parameter_meta {
        reference_fasta: "Reference genome in FASTA format"
        reference_fasta_index: "Index file for the reference genome FASTA"
        tumor_bam: "Input BAM file with aligned reads from tumor sample"
        tumor_bam_index: "Index file for the tumor BAM file"
        normal_bam: "Input BAM file with aligned reads from normal sample"
        normal_bam_index: "Index file for the normal BAM file"
        calling_regions_bed: "Optional BED file specifying regions to call variants"
        output_dir: "Directory to store Manta output"
        threads: "Number of threads to use"
        modify_disk_size_gb: "Additional disk size in GB to allocate"
    }

    input {
        File reference_fasta
        File reference_fasta_index
        File tumor_bam
        File tumor_bam_index
        File normal_bam
        File normal_bam_index
        File? calling_regions_bed
        String output_dir = "manta_output"
        Int threads = 4
        Int modify_disk_size_gb = 0
    }

    Int disk_size_gb = ceil(size(reference_fasta, "GB") * 2)
        + ceil(size(tumor_bam, "GB"))
        + ceil(size(normal_bam, "GB"))
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

        configManta.py \
            --normalBam "~{normal}" \
            --tumorBam "~{tumor}" \
            --referenceFasta "$ref_fasta" \
            ~{if defined(calling_regions_bed) then "--callRegions '" + calling_regions_bed + "'" else ""} \
            --runDir "~{output_dir}"

        "~{output_dir}/runWorkflow.py" -j "~{threads}"

        rm -rf "$ref_fasta" "$ref_fasta.fai" "~{tumor}" "~{tumor}.bai" "~{normal}" "~{normal}.bai"
    >>>

    output {
        Directory manta_output = output_dir
        File indel_candidates = "~{output_dir}/results/variants/candidateSmallIndels.vcf.gz"
        File indel_candidates_index = "~{output_dir}/results/variants/candidateSmallIndels.vcf.gz.tbi"
        File log_file = "~{output_dir}/workspace/pyflow.data/logs/pyflow_log.txt"
    }

    requirements {
        container: "quay.io/biocontainers/manta:1.6.0--py27h9948957_6"
        cpu: threads
        memory: "25 GB"
        disks: "~{disk_size_gb} GB"
    }
}
