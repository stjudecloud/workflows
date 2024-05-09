## [Homepage](https://github.com/lh3/bwa)

version 1.1

# TODO there are probably BWA params we can expose. Have not checked

task bwa_aln {
    meta {
        description: "Maps Single-End FASTQ files to BAM format using bwa aln"
        outputs: {
            bam: "Aligned BAM format file"
        }
    }

    parameter_meta {
        fastq: "Input FASTQ file to align with bwa"  # TODO verify can be gzipped or compressed
        bwa_db_tar_gz: "Gzipped tar archive of the bwa reference files. Files should be at the root of the archive."
        prefix: "Prefix for the BAM file. The extension `.bam` will be added."
        read_group: {
            description: "Read group information for BWA to insert into the header. BWA format: '@RG\tID:foo\tSM:bar'",
            common: true
        }
        use_all_cores: {
            description: "Use all cores? Recommended for cloud environments.",
            common: true
        }
        ncpu: {
            description: "Number of cores to allocate for task",
            common: true
        }
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File fastq
        File bwa_db_tar_gz
        String prefix = sub(
            basename(fastq),
            "([_\\.][rR][12])?(\\.subsampled)?\\.(fastq|fq)(\\.gz)?$",
            ""
        )
        String read_group = ""
        Boolean use_all_cores = false
        Int ncpu = 2
        Int modify_disk_size_gb = 0
    }

    String output_bam = prefix + ".bam"

    Float input_fastq_size = size(fastq, "GiB")
    Float reference_size = size(bwa_db_tar_gz, "GiB")
    Int disk_size_gb = (
        ceil((input_fastq_size + reference_size) * 2) + 10 + modify_disk_size_gb
    )

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi
        # -1 because samtools uses one more core than `--threads` specifies
        let "samtools_cores = $n_cores - 1"

        mkdir bwa_db
        tar -C bwa_db -xzf ~{bwa_db_tar_gz} --no-same-owner
        PREFIX=$(basename bwa_db/*.ann ".ann")

        bwa aln -t "$n_cores" bwa_db/"$PREFIX" ~{fastq} > sai

        bwa samse \
            ~{if read_group != "" then "-r '" + read_group + "'" else ""} \
            bwa_db/"$PREFIX" \
            sai \
            ~{fastq} \
            | samtools view --threads "$samtools_cores" -hb - \
            > ~{output_bam}

        rm -r bwa_db
    >>>

    output {
        File bam = output_bam
    }

    runtime {
        cpu: ncpu
        memory: "5 GB"
        disks: "~{disk_size_gb} GB"
        container: 'ghcr.io/stjudecloud/bwa:0.7.17-0'
        maxRetries: 1
    }
}

task bwa_aln_pe {
    meta {
        description: "Maps Paired-End FASTQ files to BAM format using bwa aln"
        outputs: {
            bam: "Aligned BAM format file"
        }
    }

    parameter_meta {
        read_one_fastq_gz: {
            description: "Input gzipped FASTQ read one file to align with bwa",
            stream: false
        }  # TODO verify can be gzipped or compressed
        read_two_fastq_gz: {
            description: "Input gzipped FASTQ read two file to align with bwa",
            stream: false
        }
        bwa_db_tar_gz: "Gzipped tar archive of the bwa reference files. Files should be at the root of the archive."
        prefix: "Prefix for the BAM file. The extension `.bam` will be added."
        read_group: {
            description: "Read group information for BWA to insert into the header. BWA format: '@RG\tID:foo\tSM:bar'",
            common: true
        }
        use_all_cores: {
            description: "Use all cores? Recommended for cloud environments.",
            common: true
        }
        ncpu: {
            description: "Number of cores to allocate for task",
            common: true
        }
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File read_one_fastq_gz
        File read_two_fastq_gz
        File bwa_db_tar_gz
        String prefix = sub(
            basename(read_one_fastq_gz),
            "([_\\.][rR][12])?(\\.subsampled)?\\.(fastq|fq)(\\.gz)?$",
            ""
        )
        String read_group = ""
        Boolean use_all_cores = false
        Int ncpu = 4
        Int modify_disk_size_gb = 0
    }

    String output_bam = prefix + ".bam"

    Float input_fastq_size = (
        size(read_one_fastq_gz, "GiB") + size(read_two_fastq_gz, "GiB")
    )
    Float reference_size = size(bwa_db_tar_gz, "GiB")
    Int disk_size_gb = (
        ceil((input_fastq_size + reference_size) * 2) + 10 + modify_disk_size_gb
    )

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi
        # -1 because samtools uses one more core than `--threads` specifies
        let "samtools_cores = $n_cores - 1"

        mkdir bwa_db
        tar -C bwa_db -xzf ~{bwa_db_tar_gz} --no-same-owner
        PREFIX=$(basename bwa_db/*.ann ".ann")

        ln -sf ~{read_one_fastq_gz}
        ln -sf ~{read_two_fastq_gz}

        bwa sampe \
            ~{if read_group != "" then "-r '"+read_group+"'" else ""} \
            bwa_db/"$PREFIX" \
            <(bwa aln -t "$n_cores" bwa_db/"$PREFIX" ~{basename(read_one_fastq_gz)}) \
            <(bwa aln -t "$n_cores" bwa_db/"$PREFIX" ~{basename(read_two_fastq_gz)}) \
            ~{basename(read_one_fastq_gz)} ~{basename(read_two_fastq_gz)} \
            | samtools view --threads "$samtools_cores" -hb - \
            > ~{output_bam}

        rm -r bwa_db
    >>>

    output {
        File bam = output_bam
    }

    runtime {
        cpu: ncpu
        memory: "17 GB"
        disks: "~{disk_size_gb} GB"
        container: 'ghcr.io/stjudecloud/bwa:0.7.17-0'
        maxRetries: 1
    }
}

task bwa_mem {
    meta {
        description: "Maps FASTQ files to BAM format using bwa mem"
        outputs: {
            bam: "Aligned BAM format file"
        }
    }

    parameter_meta {
        read_one_fastq_gz: "Input gzipped FASTQ read one file to align with bwa"  # TODO verify can be gzipped or compressed
        bwa_db_tar_gz: "Gzipped tar archive of the bwa reference files. Files should be at the root of the archive."
        read_two_fastq_gz: "Input gzipped FASTQ read two file to align with bwa"
        prefix: "Prefix for the BAM file. The extension `.bam` will be added."
        read_group: {
            description: "Read group information for BWA to insert into the header. BWA format: '@RG\tID:foo\tSM:bar'",
            common: true
        }
        use_all_cores: {
            description: "Use all cores? Recommended for cloud environments.",
            common: true
        }
        ncpu: {
            description: "Number of cores to allocate for task",
            common: true
        }
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File read_one_fastq_gz
        File bwa_db_tar_gz
        File? read_two_fastq_gz
        String prefix = sub(
            basename(read_one_fastq_gz),
            "([_\\.][rR][12])?(\\.subsampled)?\\.(fastq|fq)(\\.gz)?$",
            ""
        )
        String read_group = ""
        Boolean use_all_cores = false
        Int ncpu = 4
        Int modify_disk_size_gb = 0
    }

    String output_bam = prefix + ".bam"

    Float input_fastq_size = size(read_one_fastq_gz, "GiB") + size(read_two_fastq_gz, "GiB")
    Float reference_size = size(bwa_db_tar_gz, "GiB")
    Int disk_size_gb = (
        ceil((input_fastq_size + reference_size) * 2) + 10 + modify_disk_size_gb
    )

    File read_two_file = select_first([read_two_fastq_gz, ''])

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi
        # -1 because samtools uses one more core than `--threads` specifies
        let "samtools_cores = $n_cores - 1"

        mkdir bwa_db
        tar -C bwa_db -xzf ~{bwa_db_tar_gz} --no-same-owner
        PREFIX=$(basename bwa_db/*.ann ".ann")

        ln -sf ~{read_one_fastq_gz}
        ~{if defined(read_two_fastq_gz) then "ln -sf "+read_two_fastq_gz+"" else ""}

        bwa mem \
            -t "$n_cores" \
            ~{if read_group != "" then "-R '"+read_group+"'" else ""} \
            bwa_db/"$PREFIX" \
            ~{basename(read_one_fastq_gz)} \
            ~{basename(read_two_file)} \
            | samtools view --no-PG --threads "$samtools_cores" -hb - \
            > ~{output_bam}

        rm -r bwa_db
    >>>

    output {
        File bam = output_bam
    }

    runtime {
        cpu: ncpu
        memory: "25 GB"
        disks: "~{disk_size_gb} GB"
        container: 'ghcr.io/stjudecloud/bwa:0.7.17-0'
        maxRetries: 1
    }
}

task build_bwa_db {
    meta {
        description: "Creates a BWA index and returns it as a compressed tar archive"
        outputs: {
            bwa_db_tar_gz: "Tarballed bwa reference files"
        }
    }

    parameter_meta {
        reference_fasta: "Input reference Fasta file to index with bwa. Should be compressed with gzip."
        db_name: {
            description: "Name of the output gzipped tar archive of the bwa reference files. The extension `.tar.gz` will be added.",
            common: true
        }
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File reference_fasta
        String db_name = "bwa_db"
        Int modify_disk_size_gb = 0
    }

    Float input_fasta_size = size(reference_fasta, "GiB")
    Int disk_size_gb = ceil(input_fasta_size * 2) + 10 + modify_disk_size_gb
    String bwa_db_out_name = db_name + ".tar.gz"

    command <<<
        set -euo pipefail

        ref_fasta=~{basename(reference_fasta, ".gz")}
        gunzip -c ~{reference_fasta} > "$ref_fasta" \
            || ln -sf ~{reference_fasta} "$ref_fasta"

        bwa index "$ref_fasta"

        tar -czf ~{bwa_db_out_name} "$ref_fasta"*
    >>>

    output {
        File bwa_db_tar_gz = bwa_db_out_name
    }

    runtime {
        memory: "5 GB"
        disks: "~{disk_size_gb} GB"
        container: 'ghcr.io/stjudecloud/bwa:0.7.17-0'
        maxRetries: 1
    }
}
