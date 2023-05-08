## # BWA
##
## This WDL tool wraps [BWA](https://github.com/lh3/bwa).
## BWA aligns sequencing reads to a reference genome.

version 1.0

task bwa_aln {
    input {
        File fastq
        String output_bam = basename(fastq, ".fq.gz") + ".bam"
        File bwadb_tar_gz
        String read_group = ""
        Int ncpu = 1
        Int memory_gb = 5
        Int? disk_size_gb
        Int max_retries = 1
        Boolean use_all_cores = false
    }

    Float input_fastq_size = size(fastq, "GiB")
    Float reference_size = size(bwadb_tar_gz, "GiB")
    Int disk_size = select_first([disk_size_gb, ceil((input_fastq_size * 2) + (reference_size * 2))])

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if [ "~{use_all_cores}" = "true" ]; then
            n_cores=$(grep -c ^processor /proc/cpuinfo)
        fi

        mkdir bwa

        tar -C bwa -xzf ~{bwadb_tar_gz}
        PREFIX=$(basename bwa/*.ann ".ann")

        bwa aln -t "${n_cores}" bwa/"$PREFIX" ~{fastq} > sai

        bwa samse \
            ~{if read_group != "" then "-r '" else ""}~{read_group}~{if read_group != "" then "'" else ""} \
            bwa/"$PREFIX" sai ~{fastq} | samtools view -@ "${n_cores}" -hb - > ~{output_bam}
    >>>

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        cpu: ncpu
        docker: 'ghcr.io/stjudecloud/bwa:1.0.2'
        maxRetries: max_retries
    }

    output {
        File bam = output_bam
    }

    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
        description: "This WDL tool maps single-end fastq files to BAM format using bwa aln."
    }

    parameter_meta {
        fastq: "Input FastQ file to align with bwa"
        bwadb_tar_gz: "Gzipped tar archive of the bwa reference files. Files should be at the root of the archive."
    }
}

task bwa_aln_pe {
    input {
        File fastq1
        File fastq2
        String output_bam = basename(fastq1, ".fq.gz") + ".bam"
        File bwadb_tar_gz
        String read_group = ""
        Int ncpu = 1
        Int memory_gb = 5
        Int? disk_size_gb
        Int max_retries = 1
        Boolean use_all_cores = false
    }

    Float input_fastq_size = size(fastq1, "GiB") + size(fastq2, "GiB")
    Float reference_size = size(bwadb_tar_gz, "GiB")
    Int disk_size = select_first([disk_size_gb, ceil((input_fastq_size * 2) + (reference_size * 2))])

    command <<<
        set -xeuo pipefail

        n_cores=~{ncpu}
        if [ "~{use_all_cores}" = "true" ]; then
            n_cores=$(grep -c ^processor /proc/cpuinfo)
        fi

        mkdir bwa

        tar -C bwa -xzf ~{bwadb_tar_gz}
        PREFIX=$(basename bwa/*.ann ".ann")

        bwa aln -t "${n_cores}" bwa/"$PREFIX" ~{fastq1} > sai_1
        bwa aln -t "${n_cores}" bwa/"$PREFIX" ~{fastq2} > sai_2

        bwa sampe \
             ~{if read_group != "" then "-r '" else ""}~{read_group}~{if read_group != "" then "'" else ""} \
            bwa/"$PREFIX" sai_1 sai_2 ~{fastq1} ~{fastq2} | samtools view -@ "${n_cores}" -hb - > ~{output_bam}
    >>>

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        cpu: ncpu
        docker: 'ghcr.io/stjudecloud/bwa:1.0.2'
        maxRetries: max_retries
    }

    output {
        File bam = output_bam
    }

    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
        description: "This WDL tool maps paired-end fastq files to BAM format using bwa aln."
    }

    parameter_meta {
        fastq1: "Input FastQ read 1 file to align with bwa"
        fastq2: "Input FastQ read 2 file to align with bwa"
        bwadb_tar_gz: "Gzipped tar archive of the bwa reference files. Files should be at the root of the archive."
    }
}

task bwa_mem {
    input {
        File fastq
        String output_bam = basename(fastq, ".fq.gz") + ".bam"
        File bwadb_tar_gz
        String read_group = ""
        Int ncpu = 1
        Int memory_gb = 5
        Int? disk_size_gb
        Int max_retries = 1
        Boolean use_all_cores = false
    }

    Float input_fastq_size = size(fastq, "GiB")
    Float reference_size = size(bwadb_tar_gz, "GiB")
    Int disk_size = select_first([disk_size_gb, ceil((input_fastq_size * 2) + reference_size)])

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if [ "~{use_all_cores}" = "true" ]; then
            n_cores=$(grep -c ^processor /proc/cpuinfo)
        fi

        mkdir bwa

        tar -C bwa -xzf ~{bwadb_tar_gz}
        PREFIX=$(basename bwa/*.ann ".ann")

        bwa mem \
            -t "$n_cores" \
            ~{if read_group != "" then "-r '" else ""}~{read_group}~{if read_group != "" then "'" else ""} \
            bwa/"$PREFIX" ~{fastq} | samtools view -b - > ~{output_bam}
    >>>

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        cpu: ncpu
        docker: 'ghcr.io/stjudecloud/bwa:1.0.2'
        maxRetries: max_retries
    }

    output {
        File bam = output_bam
    }

    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
        description: "This WDL tool maps fastq files to BAM format using bwa mem."
    }

    parameter_meta {
        fastq: "Input FastQ file to align with bwa"
        bwadb_tar_gz: "Gzipped tar archive of the bwa reference files. Files should be at the root of the archive."
    }
}

task build_bwa_db {
    input {
        File reference_fasta
        String bwadb_out_name = "bwadb.tar.gz"
        Int memory_gb = 5
        Int? disk_size_gb
        Int max_retries = 1
    }

    Float input_fasta_size = size(reference_fasta, "GiB")
    Int disk_size = select_first([disk_size_gb, ceil((input_fasta_size * 2))])

    command <<<
        set -euo pipefail

        orig_fasta=~{reference_fasta}
        ref_fasta=$(basename "${orig_fasta%.gz}")
        gunzip -c ~{reference_fasta} > "$ref_fasta" || cp ~{reference_fasta} "$ref_fasta"

        bwa index "$ref_fasta"

        tar -czf ~{bwadb_out_name} "${ref_fasta}*"
    >>>

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        docker: 'ghcr.io/stjudecloud/bwa:1.0.2'
        maxRetries: max_retries
    }

    output {
        File bwadb_tar_gz = bwadb_out_name
    }

    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
        description: "This WDL tool creates a BWA index and returns it as a compressed tar archive."
    }

    parameter_meta {
        reference_fasta: "Input reference Fasta file to index with bwa. Should be compressed with gzip."
        bwadb_out_name: "Name of the output gzipped tar archive of the bwa reference files."
    }
}

task format_rg_for_bwa {
    input {
        String read_group
        Int max_retries = 1
    }

    command <<<
        echo "@RG\t~{read_group}" | sed 's/ /\\t/g' > output.txt
    >>>

    output {
        String formatted_rg = read_string("output.txt")
    }

    runtime {
        memory: "1 GB"
        disk: "1 GB"
        docker: 'ghcr.io/stjudecloud/util:1.2.0'
        maxRetries: max_retries
    }

    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
        description: "This WDL tool converts read group records from the BAM-formatted strings to strings expected by bwa."
    }

    parameter_meta {
        read_group: "Read group string"
    }
}
