version 1.0

task build_db {
    input {
        String filename = "fastq-screen-db"
        Int? max_retries = 1
    }

    String tar_filename = filename + ".tar.gz"

    command {
        fastq_screen --get_genomes 
        mv FastQ_Screen_Genomes/ ${filename}/
        tar -czf ${tar_filename} ${filename}/
    }
 #
    runtime {
        disk: "30 GB"
        docker: "fastq_screen:1.0.0-alpha"
        maxRetries: max_retries
    }

    output {
        File db = tar_filename
    }

    meta {
        author: "Clay McLeod"
        email: "clay.mcleod@STJUDE.org"
        description: ""
    }

    parameter_meta {}
} 
