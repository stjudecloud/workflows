version 1.2

task align {
    meta {
        description: "Align DNA sequences against a large reference database using BWA-MEM2"
        outputs: {
            alignments: "The output alignment file in SAM format"
        }
    }

    parameter_meta {
        read_one_fastq_gz: "Input gzipped FASTQ read one file to align with BWA-MEM2"
        reference_index: "The BWA-MEM2 index file for the reference genome"
        read_group: "The read group string to be included in the SAM header. Format: '@RG\\tID:foo\\tSM:bar'"
        read_two_fastq_gz: "Input gzipped FASTQ read two file to align with BWA-MEM2"
        prefix: "Prefix for the BAM file. The extension `.bam` will be added."
        smart_pairing: "If true, enable smart pairing mode for paired-end reads"
        skip_mate_rescue: "If true, skip mate rescue for paired-end reads"
        threads: "Number of threads to use for alignment"
        modify_disk_size_gb: "Additional disk space to allocate (in GB)"
        seed_length: "Seed value for the BWA-MEM2 aligner"
        min_score: "Minimum score threshold for reporting alignments"
    }

    input {
        File read_one_fastq_gz
        File reference_index
        String read_group
        File? read_two_fastq_gz
        String prefix = sub(
            basename(read_one_fastq_gz),
            "([_\\.][rR][12])?(\\.subsampled)?\\.(fastq|fq)(\\.gz)?$",
            ""
        )
        Boolean smart_pairing = false
        Boolean skip_mate_rescue = false
        Int threads = 4
        Int modify_disk_size_gb = 0
        Int seed_length = 19
        Int min_score = 30
    }

    String output_name = prefix + ".bam"
    Int disk_size_gb = ceil((
            size(read_one_fastq_gz, "GiB") + size(read_two_fastq_gz, "GiB")
        ) * 2)
        + ceil(size(reference_index, "GiB"))
        + 10
        + modify_disk_size_gb

    command <<<
        set -euo pipefail

        mkdir bwa_db
        tar -C bwa_db -xzf "~{reference_index}" --no-same-owner
        PREFIX=$(basename bwa_db/*.ann ".ann")

        bwa-mem2 mem \
            -t ~{threads} \
            -R "~{read_group}" \
            -k ~{seed_length} \
            -T ~{min_score} \
            ~{if smart_pairing then "-p" else ""} \
            ~{if skip_mate_rescue then "-S" else ""} \
            bwa_db/"$PREFIX" \
            "~{read_one_fastq_gz}" \
            ~{if defined(read_two_fastq_gz) then "\"~{read_two_fastq_gz}\"" else ""} |
        samtools view -b -o "~{output_name}" -

        rm -r bwa_db
    >>>

    output {
        File alignments = output_name
    }

    requirements {
        container: "ghcr.io/stjudecloud/bwamem2:branch-minimap2-2.3-0"
        cpu: threads
        memory: "64 GB"
        disks: "~{disk_size_gb} GB"
    }
}

task index {
    meta {
        description: "Index a reference genome for alignment with minimap2"
        outputs: {
            reference_index: "The minimap2 index file for the reference genome"
        }
    }

    parameter_meta {
        reference_fasta: "The reference genome in FASTA format to be indexed"
        db_name: "The base name for the output index files"
        modify_disk_size_gb: "Additional disk space to allocate (in GB)"
    }

    input {
        File reference_fasta
        String db_name = "reference"
        Int modify_disk_size_gb = 0
    }

    Float input_fasta_size = size(reference_fasta, "GiB")
    Int disk_size_gb = ceil(input_fasta_size * 2) + 10 + modify_disk_size_gb
    String bwa_db_out_name = db_name + ".tar.gz"

    command <<<
        set -euo pipefail

        ref_fasta=~{basename(reference_fasta, ".gz")}
        gunzip -c "~{reference_fasta}" > "$ref_fasta" \
            || ln -sf "~{reference_fasta}" "$ref_fasta"

        bwa-mem2 index \
            "$ref_fasta"

        tar -czf "~{bwa_db_out_name}" "$ref_fasta"*

        rm -r "$ref_fasta"
    >>>

    output {
        File reference_index = bwa_db_out_name
    }

    requirements {
        container: "ghcr.io/stjudecloud/bwamem2:branch-minimap2-2.3-0"
        cpu: 1
        memory: "120 GB"
        disks: "~{disk_size_gb} GB"
    }
}
