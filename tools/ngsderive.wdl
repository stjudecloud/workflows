## Description:
##
## This WDL tool wraps the ngsderive package (https://github.com/claymcleod/ngsderive).
## ngsderive is a utility tool to backwards compute strandedness, readlength, instrument
## for next-generation sequencing data.

version 1.0

task infer_strand {
    input {
        File bam
        File bai
        File? gtf
        Int max_retries = 1
    }

    String out_file = basename(bam, ".bam") + ".strandedness.txt"
    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil(((bam_size) * 2) + 10)
 
    command {
        sort -k1,1 -k4,4n -k5,5n ~{gtf} | bgzip > annotation.gtf.gz
        tabix -p gff annotation.gtf.gz
        mv ~{bai} ~{bam}.bai 
        ngsderive strandedness ~{bam} -g annotation.gtf.gz > ~{out_file}
        awk 'NR > 1' ~{out_file} | cut -d$'\t' -f5 > strandedness.txt
    }

    runtime {
        disk: disk_size + " GB"
        docker: 'stjudecloud/ngsderive:1.0.0-alpha'
        maxRetries: max_retries
    }

    output {
        String strandedness = read_string("strandedness.txt")
        File strandedness_file = out_file
    }
}

task instrument {
    input {
        File bam
        Int max_retries = 1
        String wait_var = ""
    }

    String out_file = basename(bam, ".bam") + ".instrument.txt"
    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil(((bam_size) * 1.25) + 10)

    command {
        ngsderive instrument ~{bam} > ~{out_file}
    }

    runtime {
        disk: disk_size + " GB"
        docker: 'stjudecloud/ngsderive:1.0.0-alpha'
        maxRetries: max_retries
    }

    output {
        File instrument_file = out_file
    }
}

task readlen {
    input {
        File bam
        Int max_retries = 1
        String wait_var = ""
    }

    String out_file = basename(bam, ".bam") + ".readlen.txt"
    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil(((bam_size) * 1.25) + 10)

    command {
        ngsderive readlen ~{bam} > ~{out_file}
    }

    runtime {
        disk: disk_size + " GB"
        docker: 'stjudecloud/ngsderive:1.0.0-alpha'
        maxRetries: max_retries
    }

    output {
        File readlen_file = out_file
    }
}
