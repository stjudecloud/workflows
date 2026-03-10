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
        output_dir: "Directory to store Manta output"
        threads: "Number of threads to use"
        modify_disk_size_gb: "Additional disk size in GB to allocate"
    }

    input {
        File reference_fasta
        File reference_fasta_index
        File bam
        File bam_index
        String output_dir = "manta_output"
        Int threads = 4
        Int modify_disk_size_gb = 0
    }

    Int disk_size_gb = ceil(size(reference_fasta, "GB") * 2)
        + ceil(size(bam, "GB"))
        + 20
        + modify_disk_size_gb

    command <<<
        set -euo pipefail

        ref_fasta=~{basename(reference_fasta, ".gz")}
        gunzip -c "~{reference_fasta}" > "$ref_fasta" \
            || ln -sf "~{reference_fasta}" "$ref_fasta"
        ln -sf "~{reference_fasta_index}" "$ref_fasta.fai"

        configManta.py \
            --bam "~{bam}" \
            --referenceFasta "$ref_fasta" \
            --runDir "~{output_dir}"

        "~{output_dir}/runWorkflow.py" -j "~{threads}"

        rm -rf "$ref_fasta" "$ref_fasta.fai"
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
            log_file: "Log file from the Manta workflow execution",
        }
    }

    parameter_meta {
        reference_fasta: "Reference genome in FASTA format"
        tumor_bam: "Input BAM file with aligned reads from tumor sample"
        tumor_bam_index: "Index file for the tumor BAM file"
        normal_bam: "Input BAM file with aligned reads from normal sample"
        normal_bam_index: "Index file for the normal BAM file"
        output_dir: "Directory to store Manta output"
        threads: "Number of threads to use"
        modify_disk_size_gb: "Additional disk size in GB to allocate"
    }

    input {
        File reference_fasta
        File tumor_bam
        File tumor_bam_index
        File normal_bam
        File normal_bam_index
        String output_dir = "manta_output"
        Int threads = 4
        Int modify_disk_size_gb = 0
    }

    Int disk_size_gb = ceil(size(reference_fasta, "GB") * 2)
        + ceil(size(tumor_bam, "GB"))
        + ceil(size(normal_bam, "GB"))
        + 20
        + modify_disk_size_gb

    command <<<
        set -euo pipefail

        ref_fasta=~{basename(reference_fasta, ".gz")}
        gunzip -c "~{reference_fasta}" > "$ref_fasta" \
            || ln -sf "~{reference_fasta}" "$ref_fasta"

        configManta.py \
            --normalBam "~{normal_bam}" \
            --tumorBam "~{tumor_bam}" \
            --referenceFasta "$ref_fasta" \
            --runDir "~{output_dir}"

        "~{output_dir}/runWorkflow.py" -j "~{threads}"

        rm -rf "$ref_fasta"
    >>>

    output {
        Directory manta_output = output_dir
        File indel_candidates = "~{output_dir}/results/variants/candidateSmallIndels.vcf.gz"
        File log_file = "~{output_dir}/workspace/pyflow.data/logs/pyflow_log.txt"
    }

    requirements {
        container: "quay.io/biocontainers/manta:1.6.0--py27h9948957_6"
        cpu: threads
        memory: "25 GB"
        disks: "~{disk_size_gb} GB"
    }
}
