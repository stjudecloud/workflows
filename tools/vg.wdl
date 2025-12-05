version 1.2

task index {
    meta {
        description: "Index a reference genome for alignment with vg giraffe"
        outputs: {
            reference_index: "The vg giraffe index file for the reference genome"
        }
    }

    parameter_meta {
        reference_fasta: "The reference genome in FASTA format to be indexed"
        vcf_files: "VCF(s) containing variants to augment the graph"
        transcript_gff: "GFF(s) containing transcript annotations"
        db_prefix: "The base name for the output index files"
        gff_feature: "The feature type in the GFF to use for transcripts"
        gff_id_tag: "The attribute tag in the GFF to use as transcript ID"
        workflow: {
            description: "The vg autoindex workflow to use",
            choices: [
                "map",
                "mpmap",
                "rpvg",
                "giraffe",
                "sr-giraffe",
                "lr-giraffe",
            ],
        }
        modify_disk_size_gb: "Additional disk space to allocate (in GB)"
        threads: "Number of threads to use for indexing"
    }

    input {
        File reference_fasta
        Array[File] vcf_files = []
        Array[File] transcript_gff = []
        String db_prefix = "reference"
        String gff_feature = "exon"
        String gff_id_tag = "transcript_id"
        String workflow = "giraffe"
        Int modify_disk_size_gb = 0
        Int threads = 4
    }

    Float input_fasta_size = size(reference_fasta, "GiB")
    Float vcf_size = size(vcf_files, "GiB")
    Float transcript_gff_size = size(transcript_gff, "GiB")
    Int disk_size_gb = ceil(input_fasta_size * 2)
        + ceil(vcf_size * 2)
        + ceil(transcript_gff_size * 2)
        + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        ref_fasta=~{basename(reference_fasta, ".gz")}
        gunzip -c "~{reference_fasta}" > "$ref_fasta" \
            || ln -sf "~{reference_fasta}" "$ref_fasta"

        vg autoindex \
            --workflow "~{workflow}" \
            -r "$ref_fasta" \
            -p "~{db_prefix}" \
            ~{sep(" ", prefix("-v ", quote(vcf_files)))} \
            ~{sep(" ", prefix("-x ", quote(transcript_gff)))} \
            -t ~{threads} \
            --gff-feature "~{gff_feature}" \
            --gff-tx-tag "~{gff_id_tag}"
    >>>

    output {
        Array[File] reference_index = glob("~{db_prefix}*")
    }

    requirements {
        container: "quay.io/biocontainers/vg:1.70.0--h9ee0642_0"
        cpu: threads
        memory: "120 GB"
        disks: "~{disk_size_gb} GB"
    }
}
