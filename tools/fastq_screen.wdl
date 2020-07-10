## # FastQ Screen
##
## Methods for bootstrapping and running [FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/_build/html/index.html)

version 1.0

task build_db {
    input {
        String filename = "fastq-screen-db"
        Int max_retries = 1
    }

    String tar_filename = filename + ".tar.gz"

    command {
        set -euo pipefail
        
        fastq_screen --get_genomes
        mv FastQ_Screen_Genomes/ ~{filename}/
        tar -czf ~{tar_filename} ~{filename}/
    }
 
    runtime {
        disk: "30 GB"
        docker: "stjudecloud/fastq_screen:1.0.0"
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
        String format
        Int num_reads = 100000
        Int max_retries = 1
    }

    Float db_size = size(db, "GiB")
    Float read1_size = size(read1, "GiB")
    Float read2_size = size(read2, "GiB")
    Int disk_size = ceil((db_size * 2) + read1_size + read2_size + 5)

    String sample_basename = basename(read1, "_R1.fastq")
    String db_name = basename(db)

    command {
        set -euo pipefail
        
        cp ~{db} /tmp
        tar -xzf /tmp/~{db_name} -C /tmp/

        cat ~{read1} ~{read2} > ~{sample_basename}.fastq

        format_arg=''
        if [[ "~{format}" = "illumina1.3" ]]; then
            format_arg='--illumina1_3'
        fi;
        fastq_screen \
            $format_arg \
            --subset ~{num_reads} \
            --conf /home/fastq_screen.conf \
            --aligner bowtie2 \
            ~{sample_basename}.fastq
    }
 
    runtime {
        disk: disk_size + " GB"
        docker: "stjudecloud/fastq_screen:1.0.0"
        maxRetries: max_retries
    }

    output {
        Array[File] out_files = glob("~{sample_basename}_screen*")
    }

    meta {
        author: "Andrew Frantz"
        email: "andrew.frantz@stjude.org"
        description: "This WDL tool runs FastQ Screen on a pair of fastqs."
    }
}
