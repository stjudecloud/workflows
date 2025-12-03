version 1.2

task align {
    meta {
        description: "Align DNA or mRNA sequences against a large reference database"
        outputs: {
            alignments: "The output alignment file in SAM or PAF format"
        }
    }

    parameter_meta {
        reads: "The input reads file in FASTQ format"
        reference_index: "The minimap2 index file for the reference genome"
        read_group: "The read group string to be included in the SAM header. Format: '@RG\\tID:foo\\tSM:bar'"
        output_name: "The name of the output alignment file"
        output_paf: "If true, output in PAF format instead of SAM"
        cigar_in_paf: "If true and outputting PAF, include CIGAR strings in the PAF output"
        ignore_base_quality: "If true, ignore base quality scores during alignment"
        output_md_tag: "If true, include MD tags in the SAM output"
        eqx: "If true, use =/X CIGAR operators instead of M"
        soft_clip: "If true, use soft clipping for secondary alignments in SAM format"
        secondary_alignments: "If true, report secondary alignments"
        seed: "Seed value for the minimap2 aligner"
        threads: "Number of threads to use for alignment"
    }

    input {
        File reads
        File reference_index
        String read_group
        String output_name = "aligned.sam"
        Boolean output_paf = false
        Boolean cigar_in_paf = true
        Boolean ignore_base_quality = false
        Boolean output_md_tag = true
        Boolean eqx = false
        Boolean soft_clip = true
        Boolean secondary_alignments = true
        Int seed = 11
        Int threads = 3
    }

    command <<<
        minimap2 \
            ~{if output_paf then "" else "-a"} \
            ~{if output_paf && cigar_in_paf then "-c" else ""} \
            ~{if ignore_base_quality then "-Q" else ""} \
            ~{if output_md_tag then "--MD" else ""} \
            ~{if eqx then "-X" else ""} \
            ~{if soft_clip then "-Y" else ""} \
            ~{if secondary_alignments then "--secondary=yes" else "--secondary=no"} \
            -t ~{threads} \
            --seed ~{seed} \
            -R "~{read_group}" \
            "~{reference_index}" \
            "~{reads}" \
            > "~{output_name}"
    >>>

    output {
        File alignments = output_name
    }

    requirements {
        container: "quay.io/biocontainers/minimap2:2.30--h577a1d6_0"
        cpu: threads
        memory: "4 GB"
    }
}

task index {
    meta {
        description: "Create a minimap2 index for a reference genome"
        outputs: {
            reference_index: "The generated minimap2 index file"
        }
    }

    parameter_meta {
        reference_fasta: "The reference genome in FASTA format to be indexed"
        alt_contigs: "Optional file containing a list of alternative contigs"
        index_name: "The name of the output index file"
        minimizer_kmer_size: "K-mer size for minimizer indexing"
        minimizer_window_size: "Window size for minimizer indexing"
    }

    input {
        File reference_fasta
        File? alt_contigs
        String index_name = "reference.mmi"
        Int minimizer_kmer_size = 15
        Int minimizer_window_size = 10
    }

    command <<<
        set -euo pipefail

        ref_fasta=~{basename(reference_fasta, ".gz")}
        gunzip -c "~{reference_fasta}" > "$ref_fasta" \
            || ln -sf "~{reference_fasta}" "$ref_fasta"

        minimap2 \
            -k ~{minimizer_kmer_size} \
            -w ~{minimizer_window_size} \
            ~{if defined(alt_contigs) then "--alt \"~{alt_contigs}\"" else ""} \
            -d "~{index_name}" \
            "$ref_fasta"
    >>>

    output {
        File reference_index = index_name
    }

    requirements {
        container: "quay.io/biocontainers/minimap2:2.30--h577a1d6_0"
        cpu: 1
        memory: "4 GB"
    }
}
