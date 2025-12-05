version 1.2

task giraffe {
    meta {
        description: "Align DNA sequences against a variation graph using vg giraffe"
        outputs: {
            alignments: "The output alignment file in GAM format"
        }
    }

    parameter_meta {
        read_one_fastq_gz: "Input gzipped FASTQ read one file to align with vg giraffe"
        gbz_graph: "The vg GBZ graph file for the reference genome"
        minimizer_index: "The vg minimizer index file for the reference genome"
        zipcode_name: "The vg zipcode name file for the reference genome"
        distance_index: "The vg distance index file for the reference genome"
        read_two_fastq_gz: "Input gzipped FASTQ read two file to align with vg giraffe"
        haploytype: "The haplotype information file"
        kff: "The KFF file containing kmer counts"
        sample_name: "The sample name to include"
        read_group: "The read group"
        output_name: "The name of the output alignment file"
        output_format: {
            description: "The output format for alignments",
            options: [
                "gam",
                "gaf",
                "json",
                "tsv",
                "SAM",
                "BAM",
                "CRAM",
            ],
        }
        preset: {
            description: "vg giraffe preset for alignment",
            options: [
                "chaining-sr",
                "default",
                "fast",
                "hifi",
                "r10",
                "srold",
            ],
        }
        threads: "Number of threads to use for alignment"
        modify_disk_size_gb: "Additional disk space to allocate (in GB)"
    }

    input {
        File read_one_fastq_gz
        File gbz_graph
        File minimizer_index
        File zipcode_name
        File distance_index
        File? read_two_fastq_gz
        File? haploytype
        File? kff
        String? sample_name
        String? read_group
        String output_name = "aligned.bam"
        String output_format = "BAM"
        String preset = "default"
        Int threads = 4
        Int modify_disk_size_gb = 0
    }

    Int disk_size_gb = ceil((
            size(read_one_fastq_gz, "GiB") + size(read_two_fastq_gz, "GiB")
        ) * 2)
        + ceil(size(gbz_graph, "GiB"))
        + ceil(size(minimizer_index, "GiB"))
        + ceil(size(distance_index, "GiB"))
        + ceil(size(zipcode_name, "GiB"))
        + 10
        + modify_disk_size_gb

    command <<<
        vg giraffe \
            -t ~{threads} \
            -Z "~{gbz_graph}" \
            -m "~{minimizer_index}" \
            -d "~{distance_index}" \
            -z "~{zipcode_name}" \
            -f "~{read_one_fastq_gz}" \
            ~{if defined(read_two_fastq_gz) then "-f \"~{read_two_fastq_gz}\"" else ""} \
            -o "~{output_format}" \
            ~{if defined(sample_name) then "--sample \"~{sample_name}\"" else ""} \
            ~{if defined(read_group) then "--read-group \"~{read_group}\"" else ""} \
            ~{if defined(haploytype) then "--haplotype-name \"~{haploytype}\"" else ""} \
            ~{if defined(kff) then "--kff-name \"~{kff}\"" else ""} \
            --parameter-preset "~{preset}" \
            > "~{output_name}"
    >>>

    output {
        File alignments = "~{output_name}"
    }

    requirements {
        container: "quay.io/biocontainers/vg:1.70.0--h9ee0642_0"
        cpu: threads
        memory: "120 GB"
        disks: "~{disk_size_gb} GB"
    }
}

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
