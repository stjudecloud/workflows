version 1.0

task build_db {
    input {
        String filename = "fastq-screen-db"
        Int max_retries = 1
    }

    String tar_filename = filename + ".tar.gz"

    command {
        fastq_screen --get_genomes
        mv FastQ_Screen_Genomes/ ~{filename}/
        tar -czf ~{tar_filename} ~{filename}/
    }
 
    runtime {
        disk: "30 GB"
        docker: "stjudecloud/fastq_screen:1.0.0-alpha"
        maxRetries: max_retries
    }

    output {
        File db = tar_filename
    }

    meta {
        author: "Clay McLeod, Andrew Frantz"
        email: "clay.mcleod@STJUDE.org, andrew.frantz@stjude.org"
        description: "This WDL tool downloads the Fast Q Screen database and archives it."
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

    String output_basename = basename(read1, "R1.fastq") + "screen"
    String db_name = basename(db)

    command {
        cp ~{db} /tmp
        tar -xsf /tmp/~{db_name} -C /tmp/

        format_arg=''
        if [[ "~{format}" = "illumina1.3" ]]; then
            format_arg='--illumina1_3'
        fi;
        fastq_screen \
            $format_arg \
            --subset ~{num_reads} \
            --conf /home/fastq_screen.conf \
            --aligner bowtie2 \
            ~{read1} ~{read2}
    }
 
    runtime {
        docker: "stjudecloud/fastq_screen:1.0.0-alpha"
        maxRetries: max_retries
    }

    output {
        Array[File] out_files = glob("~{output_basename}*")
    }

    meta {
        author: "Andrew Frantz"
        email: "andrew.frantz@stjude.org"
        description: "This WDL tool runs Fast Q Screen."
    }
}
