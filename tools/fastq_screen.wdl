## # FastQ Screen
##
## Methods for bootstrapping and running [FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/_build/html/index.html)

version 1.0

task build_db {
    input {
        String tar_filename = "fastq-screen-db.tar.gz"
        Int max_retries = 1
    }

    command {
        set -euo pipefail
        
        fastq_screen --get_genomes
        mv FastQ_Screen_Genomes/fastq_screen.conf FastQ_Screen_Genomes/fastq_screen.conf.template
        tar -czf ~{tar_filename} -C FastQ_Screen_Genomes/ .
    }
 
    runtime {
        disk: "30 GB"
        docker: 'ghcr.io/stjudecloud/fastq_screen:1.1.3'
        maxRetries: max_retries
    }

    output {
        File db = tar_filename
    }

    meta {
        author: "Clay McLeod, Andrew Frantz"
        email: "clay.mcleod@STJUDE.org, andrew.frantz@stjude.org"
        description: "This WDL tool downloads the FastQ Screen database and archives it."
    }
} 

task fastq_screen {
    input {
        File read1
        File read2
        File db
        String provided_encoding
        String inferred_encoding = ""
        Int num_reads = 0
        String? sample_name
        Int max_retries = 1
    }

    parameter_meta {
        db: "Database for FastQ Screen. Must untar directly to the genome directories."
    }

    Float db_size = size(db, "GiB")
    Float read1_size = size(read1, "GiB")
    Float read2_size = size(read2, "GiB")
    Int disk_size = ceil((db_size * 2) + read1_size + read2_size + 5)

    String inferred_basename = basename(read1, "_R1.fastq.gz")
    String sample_basename = select_first([sample_name, inferred_basename])
    String out_directory = sample_basename + ".screen"
    String out_tar_gz = out_directory + ".tar.gz"

    String format_arg = if (provided_encoding != "") then
                            if (provided_encoding == "illumina1.3") then "--illumina1_3"
                            else ""
                        else
                            if (inferred_encoding == "Illumina 1.3") then "--illumina1_3"
                            else ""

    command <<<
        set -euo pipefail

        export LC_ALL=C.UTF-8
        export LANG=C.UTF-8

        mkdir -p /tmp/FastQ_Screen_Genomes/
        tar -xzf ~{db} -C /tmp/FastQ_Screen_Genomes/ --no-same-owner

        gunzip -c ~{read1} ~{read2} > ~{sample_basename}.fastq

        fastq_screen \
            ~{format_arg} \
            --subset ~{num_reads} \
            --conf /home/fastq_screen.conf \
            --aligner bowtie2 \
            ~{sample_basename}.fastq 2>&1 \
            | sed '/Skipping DATABASE/q1;/ERR/q1' 1>&2 \
            || exit 42

        mkdir ~{out_directory}
        mv "~{out_directory}".* "~{out_directory}"
        tar -czf ~{out_tar_gz} ~{out_directory}
    >>>
 
    runtime {
        memory: "10 GB"
        disk: disk_size + " GB"
        docker: 'ghcr.io/stjudecloud/fastq_screen:1.1.3'
        maxRetries: max_retries
    }

    output {
        File results = out_tar_gz
    }

    meta {
        author: "Andrew Frantz"
        email: "andrew.frantz@stjude.org"
        description: "This WDL tool runs FastQ Screen on a sample. Exit code 42 indicates a rare intermittent bug. Job should succeed upon resubmission."
    }
}
