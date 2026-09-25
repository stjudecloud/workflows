version 1.3

enum OutputFormat {
    gam,
    gaf,
    json,
    tsv,
    SAM,
    BAM,
    CRAM,
}

# [vg giraffe presets source](https://github.com/vgteam/vg/blob/v1.70.0/src/subcommand/giraffe_main.cpp)
#
# | `-b`          | Use case                                               | Notes                                                                                                   |
# | ------------- | ------------------------------------------------------ | ------------------------------------------------------------------------------------------------------- |
# | `default`     | General-purpose short/long read mapping                | Balanced default parameters.                                                                            |
# | `fast`        | Faster mapping with reduced sensitivity                | Lowers hit caps and multimapping limits for speed.                                                      |
# | `hifi`        | PacBio HiFi long reads                                 | Uses chaining-based alignment (`align-from-chains`) tuned for HiFi error profiles.                      |
# | `r10`         | Oxford Nanopore R10 long reads                         | Uses chaining-based alignment (`align-from-chains`) tuned for R10 error profiles.                       |
# | `chaining-sr` | Short reads using the chaining-based codepath          | Shares the chaining algorithm introduced for long-read presets (`hifi`/`r10`), adapted for short reads. |
# | `srold`       | Short reads using the original (pre-chaining) codepath | Legacy short-read algorithm predating `chaining-sr`.                                                    |
enum GiraffePreset[String] {
    chaining_sr = "chaining-sr",
    default,
    fast,
    hifi,
    r10,
    srold,
}

# [vg autoindex source](https://github.com/vgteam/vg/blob/v1.70.0/src/index_registry.hpp#L139-L148)
#
# | `-w`                     | Indexes produced                                             | Notes                                         |
# | ------------------------ | ------------------------------------------------------------ | --------------------------------------------- |
# | `map`                    | Indexes for `vg map`                                         | The default workflow.                         |
# | `mpmap`                  | Indexes for `vg mpmap` (multipath mapper)                    |                                               |
# | `rpvg`                   | Indexes for `rpvg` haplotype-based transcript quantification |                                               |
# | `giraffe` / `sr-giraffe` | Indexes for `vg giraffe` on short reads                      | `giraffe` is a legacy alias for `sr-giraffe`. |
# | `lr-giraffe`             | Indexes for `vg giraffe` on long reads                       |                                               |
enum AutoindexWorkflow[String] {
    map,
    mpmap,
    rpvg,
    giraffe,
    sr_giraffe = "sr-giraffe",
    lr_giraffe = "lr-giraffe",
}

task giraffe {
    meta {
        description: "Align DNA sequences against a variation graph using vg giraffe"
        outputs: {
            alignments: "The output alignment file in GAM format",
        }
    }

    parameter_meta {
        read_one_fastq_gz: "Input gzipped FASTQ read one file to align with vg giraffe"
        gbz_graph: "The vg GBZ graph file for the reference genome"
        minimizer_index: "The vg minimizer index file for the reference genome"
        zipcode_name: "The vg zipcode name file for the reference genome"
        distance_index: "The vg distance index file for the reference genome"
        read_two_fastq_gz: "Input gzipped FASTQ read two file to align with vg giraffe"
        haplotype: "The haplotype information file"
        kff: "The KFF file containing kmer counts"
        sample_name: "The sample name to include"
        read_group: "The read group"
        output_name: "The name of the output alignment file"
        output_format: "The output format for alignments"
        preset: "vg giraffe preset for alignment"
        ncpu: "Number of threads to use for alignment"
        modify_disk_size_gb: "Additional disk space to allocate (in GB)"
    }

    input {
        File read_one_fastq_gz
        File gbz_graph
        File minimizer_index
        File zipcode_name
        File distance_index
        File? read_two_fastq_gz
        File? haplotype
        File? kff
        String? sample_name
        String? read_group
        GiraffePreset preset = GiraffePreset.default
        OutputFormat output_format = OutputFormat.BAM
        String output_name = "aligned.bam"
        Int ncpu = 4
        Int modify_disk_size_gb = 0
    }

    Int disk_size_gb = ceil((size(read_one_fastq_gz, "GB") + size(read_two_fastq_gz, "GB")
    ) * 2) + ceil(size(gbz_graph, "GB")) + ceil(size(minimizer_index, "GB")) + ceil(size(
        distance_index, "GB"
    )) + ceil(size(zipcode_name, "GB")) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail
        vg giraffe \
            -t ~{ncpu} \
            -Z "~{gbz_graph}" \
            -m "~{minimizer_index}" \
            -d "~{distance_index}" \
            -z "~{zipcode_name}" \
            -f "~{read_one_fastq_gz}" \
            ~{if defined(read_two_fastq_gz) then "-f \"~{read_two_fastq_gz}\"" else ""} \
            -o "~{output_format}" \
            ~{if defined(sample_name) then "--sample \"~{sample_name}\"" else ""} \
            ~{if defined(read_group) then "--read-group \"~{read_group}\"" else ""} \
            ~{if defined(haplotype) then "--haplotype-name \"~{haplotype}\"" else ""} \
            ~{if defined(kff) then "--kff-name \"~{kff}\"" else ""} \
            --parameter-preset "~{preset}" \
            > "~{output_name}"
    >>>

    output {
        File alignments = "~{output_name}"
    }

    requirements {
        container: "quay.io/biocontainers/vg:1.70.0--h9ee0642_0"
        cpu: ncpu
        memory: "60 GB"
        disks: "~{disk_size_gb} GB"
    }
}

task index {
    meta {
        description: "Index a reference genome for alignment with vg giraffe"
        outputs: {
            reference_index: "The vg giraffe index file for the reference genome",
        }
    }

    parameter_meta {
        reference_fasta: "The reference genome in FASTA format to be indexed"
        vcf_files: "VCF(s) containing variants to augment the graph"
        transcript_gff: "GFF(s) containing transcript annotations"
        db_prefix: "The base name for the output index files"
        gff_feature: "The feature type in the GFF to use for transcripts"
        gff_id_tag: "The attribute tag in the GFF to use as transcript ID"
        autoindex_workflow: "The vg autoindex workflow to use"
        modify_disk_size_gb: "Additional disk space to allocate (in GB)"
        ncpu: "Number of threads to use for indexing"
    }

    input {
        File reference_fasta
        Array[File] vcf_files = []
        Array[File] transcript_gff = []
        AutoindexWorkflow autoindex_workflow = AutoindexWorkflow.giraffe
        String db_prefix = "reference"
        String gff_feature = "exon"
        String gff_id_tag = "transcript_id"
        Int modify_disk_size_gb = 0
        Int ncpu = 4
    }

    Float input_fasta_size = size(reference_fasta, "GB")
    Float vcf_size = size(vcf_files, "GB")
    Float transcript_gff_size = size(transcript_gff, "GB")
    Int disk_size_gb = ceil(input_fasta_size * 2) + ceil(vcf_size * 2) + ceil(
        transcript_gff_size * 2
    ) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        ref_fasta=~{basename(reference_fasta, ".gz")}
        gunzip -c "~{reference_fasta}" > "$ref_fasta" \
            || ln -sf "~{reference_fasta}" "$ref_fasta"

        vg autoindex \
            --workflow "~{autoindex_workflow}" \
            -r "$ref_fasta" \
            -p "~{db_prefix}" \
            ~{sep(" ", prefix("-v ", quote(vcf_files)))} \
            ~{sep(" ", prefix("-x ", quote(transcript_gff)))} \
            -t ~{ncpu} \
            --gff-feature "~{gff_feature}" \
            --gff-tx-tag "~{gff_id_tag}"
    >>>

    output {
        Array[File] reference_index = glob("~{db_prefix}*")
    }

    requirements {
        container: "quay.io/biocontainers/vg:1.70.0--h9ee0642_0"
        cpu: ncpu
        memory: "120 GB"
        disks: "~{disk_size_gb} GB"
    }
}
