version 1.0

task build_db {
    input {
        String out_filename = "fastq-screen-db.tar.gz"
    }

    command {
        fastq_screen --get_genomes
        mv Fastq_Screen_Genomes/ db/
        tar -czf ${out_filename} db/
    }

    runtime {
        disk: "30 GB"
        docker: "stjudecloud/fastq_screen:bleeding-edge"
        maxRetries: max_retries
    }

    output {
        File db = $out_filename
    }

    meta {
        author: "Clay McLeod"
        email: "clay.mcleod@STJUDE.org"
        description: ""
    }

    parameter_meta {}
} 
