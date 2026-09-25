version 1.3

# [minimap2 preset options](https://lh3.github.io/minimap2/minimap2.html#8)
#
# | `-x`        | Use case                                                  | Notes                                                                                                 |
# | ----------- | --------------------------------------------------------- | ----------------------------------------------------------------------------------------------------- |
# | `map-ont`   | Noisy long reads (~10% error) to a reference genome       | Default mode.                                                                                         |
# | `lr:hq`     | Accurate long reads (<1% error) to a reference genome     | Recommended for Nanopore v14 chemistry reads (~99% accuracy); works better than `map-hifi` for these. |
# | `map-hifi`  | PacBio HiFi reads to a reference genome                   | Same as `lr:hq` but with different scoring.                                                           |
# | `map-pb`    | Older PacBio CLR reads to a reference genome              | Effectively deprecated by HiFi; prefer `map-hifi`/`lr:hq` unless working with old data.               |
# | `map-iclr`  | Illumina Complete Long Reads (ICLR) to a reference genome | Recommended by Illumina developers.                                                                   |
# | `asm5`      | Long assembly-to-reference mapping                        | For divergence not much higher than 0.1%.                                                             |
# | `asm10`     | Long assembly-to-reference mapping                        | For divergence around 1%.                                                                             |
# | `asm20`     | Long assembly-to-reference mapping                        | For divergence around several percent.                                                                |
# | `splice`    | Long-read spliced alignment                               | Long deletions become introns (`N` CIGAR op); long insertions disabled.                               |
# | `splice:hq` | Spliced alignment for accurate long RNA-seq reads         | e.g. PacBio Iso-Seq.                                                                                  |
# | `splice:sr` | Spliced alignment for short RNA-seq reads                 |                                                                                                       |
# | `sr`        | Short-read alignment without splicing                     |                                                                                                       |
# | `ava-pb`    | PacBio CLR all-vs-all overlap mapping                     |                                                                                                       |
# | `ava-ont`   | Oxford Nanopore all-vs-all overlap mapping                |                                                                                                       |
enum Preset[String] {
    sr,
    map_ont = "map-ont",
    lr_hq = "lr:hq",
    map_hifi = "map-hifi",
    map_pb = "map-pb",
    map_iclr = "map-iclr",
    asm5,
    asm10,
    asm20,
    splice,
    splice_hq = "splice:hq",
    splice_sr = "splice:sr",
    ava_pb = "ava-pb",
    ava_ont = "ava-ont",
}

task align {
    meta {
        description: "Align DNA or mRNA sequences against a large reference database"
        outputs: {
            alignments: "The output alignment file in SAM or PAF format",
        }
    }

    parameter_meta {
        read_one_fastq_gz: "Input gzipped FASTQ read one file to align with minimap2"
        reference_index: "The minimap2 index file for the reference genome"
        read_group: "The read group string to be included in the SAM header. Format: '@RG\\tID:foo\\tSM:bar'"
        read_two_fastq_gz: "Input gzipped FASTQ read two file to align with minimap2"
        preset: {
            description: "Minimap2 preset for alignment",
            external_help: "https://lh3.github.io/minimap2/minimap2.html#8",
        }
        output_name: "The name of the output alignment file"
        output_paf: "If true, output in PAF format instead of BAM"
        cigar_in_paf: "If true and outputting PAF, include CIGAR strings in the PAF output"
        ignore_base_quality: "If true, ignore base quality scores during alignment"
        output_md_tag: "If true, include MD tags in the SAM output"
        eqx: "If true, use =/X CIGAR operators instead of M"
        soft_clip: "If true, use soft clipping for secondary alignments in SAM format"
        secondary_alignments: "If true, report secondary alignments"
        seed: "Seed value for the minimap2 aligner"
        ncpu: "Number of threads to use for alignment"
        modify_disk_size_gb: "Additional disk space to allocate (in GB)"
    }

    input {
        File read_one_fastq_gz
        File reference_index
        String read_group
        File? read_two_fastq_gz
        Preset? preset = Preset.sr
        String output_name = "aligned.bam"
        Boolean output_paf = false
        Boolean cigar_in_paf = true
        Boolean ignore_base_quality = false
        Boolean output_md_tag = true
        Boolean eqx = false
        Boolean soft_clip = false
        Boolean secondary_alignments = false
        Int seed = 11
        Int ncpu = 3
        Int modify_disk_size_gb = 0
    }

    Int disk_size_gb = ceil((size(read_one_fastq_gz, "GB") + size(read_two_fastq_gz, "GB")
    ) * 3) + ceil(size(reference_index, "GB")) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        minimap2 \
            ~{if defined(preset) then "-x \"~{preset}\"" else ""} \
            ~{if output_paf then "" else "-a"} \
            ~{if output_paf && cigar_in_paf then "-c" else ""} \
            ~{if ignore_base_quality then "-Q" else ""} \
            ~{if output_md_tag then "--MD" else ""} \
            ~{if eqx then "-X" else ""} \
            ~{if soft_clip then "-Y" else ""} \
            ~{if secondary_alignments then "--secondary=yes" else "--secondary=no"} \
            -t ~{ncpu} \
            --seed ~{seed} \
            -R "~{read_group}" \
            "~{reference_index}" \
            "~{read_one_fastq_gz}" \
            ~{if defined(read_two_fastq_gz) then "\"~{read_two_fastq_gz}\"" else ""} \
            | if ~{output_paf}; then
                cat - > "~{output_name}"
            else
                samtools view -b - > "~{output_name}"
            fi
    >>>

    output {
        File alignments = output_name
    }

    requirements {
        container: "ghcr.io/stjudecloud/minimap2:2.30-0"
        cpu: ncpu
        memory: "16 GB"
        disks: "~{disk_size_gb} GB"
    }
}

task index {
    meta {
        description: "Create a minimap2 index for a reference genome"
        outputs: {
            reference_index: "The generated minimap2 index file",
        }
    }

    parameter_meta {
        reference_fasta: "The reference genome in FASTA format to be indexed"
        alt_contigs: "Optional file containing a list of alternative contigs"
        index_name: "The name of the output index file"
        minimizer_kmer_size: "K-mer size for minimizer indexing"
        minimizer_window_size: "Window size for minimizer indexing"
        ncpu: "Number of threads to use for indexing"
        modify_disk_size_gb: "Additional disk space to allocate (in GB)"
    }

    input {
        File reference_fasta
        File? alt_contigs
        String index_name = "reference.mmi"
        Int minimizer_kmer_size = 15
        Int minimizer_window_size = 10
        Int ncpu = 3
        Int modify_disk_size_gb = 0
    }

    Int disk_size_gb = ceil(size(reference_fasta, "GB")) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        ref_fasta=~{basename(reference_fasta, ".gz")}
        gunzip -c "~{reference_fasta}" > "$ref_fasta" \
            || ln -sf "~{reference_fasta}" "$ref_fasta"

        minimap2 \
            -k ~{minimizer_kmer_size} \
            -w ~{minimizer_window_size} \
            ~{if defined(alt_contigs) then "--alt \"~{alt_contigs}\"" else ""} \
            -t ~{ncpu} \
            -d "~{index_name}" \
            "$ref_fasta"

        rm -r "$ref_fasta"
    >>>

    output {
        File reference_index = index_name
    }

    requirements {
        container: "ghcr.io/stjudecloud/minimap2:2.30-0"
        cpu: ncpu
        memory: "16 GB"
        disks: "~{disk_size_gb} GB"
    }
}
