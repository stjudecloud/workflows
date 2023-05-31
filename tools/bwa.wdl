## # BWA
##
## This WDL file wraps [BWA](https://github.com/lh3/bwa).
## BWA aligns sequencing reads to a reference genome.

version 1.0

task bwa_aln {
    meta {
        description: "This WDL task maps single-end FastQ files to BAM format using bwa aln."
    }

    parameter_meta {
        fastq: "Input FastQ file to align with bwa"
        bwa_db_tar_gz: "Gzipped tar archive of the bwa reference files. Files should be at the root of the archive."
    }

    input {
        File fastq
        File bwa_db_tar_gz
        String prefix = basename(fastq, ".fq.gz")  # TODO is this right?
        String read_group = ""
        Boolean use_all_cores = false
        Int memory_gb = 5
        Int modify_disk_size_gb = 0
        Int ncpu = 1
        Int max_retries = 1
    }

    String output_bam = prefix + ".bam"

    Float input_fastq_size = size(fastq, "GiB")
    Float reference_size = size(bwa_db_tar_gz, "GiB")
    Int disk_size_gb = ceil(
        (input_fastq_size + reference_size) * 2
    ) + modify_disk_size_gb

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi

        mkdir bwa_db
        tar -C bwa_db -xzf ~{bwa_db_tar_gz} --no-same-owner
        PREFIX=$(basename bwa_db/*.ann ".ann")

        bwa aln -t "$n_cores" bwa_db/"$PREFIX" ~{fastq} > sai

        bwa samse \
            ~{if read_group != "" then "-r '" + read_group + "'" else ""} \
            bwa_db/"$PREFIX" \
            sai \
            ~{fastq} \
            | samtools view -@ "$n_cores" -hb - \
            > ~{output_bam}

        rm -r bwa_db
    >>>

    output {
        File bam = output_bam
    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        cpu: ncpu
        docker: 'ghcr.io/stjudecloud/bwa:0.7.17-0'
        maxRetries: max_retries
    }
}

task bwa_aln_pe {
    meta {
        description: "This WDL task maps paired-end FastQ files to BAM format using bwa aln."
    }

    parameter_meta {
        read_one_fastq_gz: "Input FastQ read 1 file to align with bwa"
        read_two_fastq_gz: "Input FastQ read 2 file to align with bwa"
        bwa_db_tar_gz: "Gzipped tar archive of the bwa reference files. Files should be at the root of the archive."
    }

    input {
        File read_one_fastq_gz
        File read_two_fastq_gz
        File bwa_db_tar_gz
        String prefix = basename(read_one_fastq_gz, ".fq.gz")  # TODO is this right?
        String read_group = ""
        Boolean use_all_cores = false
        Int memory_gb = 5
        Int modify_disk_size_gb = 0
        Int ncpu = 1
        Int max_retries = 1
    }

    String output_bam = prefix + ".bam"

    Float input_fastq_size = size(read_one_fastq_gz, "GiB") + size(read_two_fastq_gz, "GiB")
    Float reference_size = size(bwa_db_tar_gz, "GiB")
    Int disk_size_gb = ceil(
        (input_fastq_size + reference_size) * 2
    ) + modify_disk_size_gb

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi

        mkdir bwa_db
        tar -C bwa_db -xzf ~{bwa_db_tar_gz} --no-same-owner
        PREFIX=$(basename bwa_db/*.ann ".ann")

        bwa aln -t "$n_cores" bwa_db/"$PREFIX" ~{read_one_fastq_gz} > sai_1
        bwa aln -t "$n_cores" bwa_db/"$PREFIX" ~{read_two_fastq_gz} > sai_2

        bwa sampe \
            ~{if read_group != "" then "-r '"+read_group+"'" else ""} \
            bwa_db/"$PREFIX" \
            sai_1 sai_2 \
            ~{read_one_fastq_gz} ~{read_two_fastq_gz} \
            | samtools view -@ "$n_cores" -hb - \
            > ~{output_bam}

        rm -r bwa_db
    >>>

    output {
        File bam = output_bam
    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        cpu: ncpu
        docker: 'ghcr.io/stjudecloud/bwa:0.7.17-0'
        maxRetries: max_retries
    }
}

task bwa_mem {
    meta {
        description: "This WDL task maps FastQ files to BAM format using bwa mem."
    }

    parameter_meta {
        fastq: "Input FastQ file to align with bwa"
        bwa_db_tar_gz: "Gzipped tar archive of the bwa reference files. Files should be at the root of the archive."
    }

    input {
        File fastq
        File bwa_db_tar_gz
        String prefix = basename(fastq, ".fq.gz")  # TODO is this right?
        String read_group = ""
        Boolean use_all_cores = false
        Int memory_gb = 5
        Int modify_disk_size_gb
        Int ncpu = 1
        Int max_retries = 1
    }

    String output_bam = prefix + ".bam"

    Float input_fastq_size = size(fastq, "GiB")
    Float reference_size = size(bwa_db_tar_gz, "GiB")
    Int disk_size_gb = ceil(
        (input_fastq_size + reference_size) * 2
    ) + modify_disk_size_gb

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi

        mkdir bwa_db
        tar -C bwa_db -xzf ~{bwa_db_tar_gz} --no-same-owner
        PREFIX=$(basename bwa_db/*.ann ".ann")

        bwa mem \
            -t "$n_cores" \
            ~{if read_group != "" then "-r '"+read_group+"'" else ""} \
            bwa_db/"$PREFIX" \
            ~{fastq} \
            | samtools view -@ "$n_cores" -hb - \
            > ~{output_bam}

        rm -r bwa_db
    >>>

    output {
        File bam = output_bam
    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        cpu: ncpu
        docker: 'ghcr.io/stjudecloud/bwa:0.7.17-0'
        maxRetries: max_retries
    }
}

task build_bwa_db {
    meta {
        description: "This WDL task creates a BWA index and returns it as a compressed tar archive."
    }

    parameter_meta {
        reference_fasta: "Input reference Fasta file to index with bwa. Should be compressed with gzip."
        db_name: "Name of the output gzipped tar archive of the bwa reference files."
    }

    input {
        File reference_fasta
        String db_name = "bwa_db"
        Int memory_gb = 5
        Int modify_disk_size_gb = 0
        Int max_retries = 1
    }

    Float input_fasta_size = size(reference_fasta, "GiB")
    Int disk_size_gb = ceil(input_fasta_size * 2) + modify_disk_size_gb
    String bwa_db_out_name = db_name + ".tar.gz"

    command <<<
        set -euo pipefail

        orig_fasta=~{reference_fasta}
        ref_fasta=$(basename "${orig_fasta%.gz}")
        gunzip -c ~{reference_fasta} > "$ref_fasta" \
            || ln -s ~{reference_fasta} "$ref_fasta"

        bwa index "$ref_fasta"

        tar -czf ~{bwa_db_out_name} "${ref_fasta}"*
    >>>

    output {
        File bwa_db_tar_gz = bwa_db_out_name
    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        docker: 'ghcr.io/stjudecloud/bwa:0.7.17-0'
        maxRetries: max_retries
    }
}

task format_rg_for_bwa {
    meta {
        description: "This WDL task converts read group records from the BAM-formatted strings to strings expected by bwa."
    }

    parameter_meta {
        read_group: "Read group string"
    }

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
}
