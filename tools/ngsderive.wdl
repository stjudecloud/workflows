## # ngsderive
##
## This WDL tool wraps the [ngsderive package](https://github.com/claymcleod/ngsderive).
## ngsderive is a utility tool to backwards compute strandedness, readlength, instrument
## for next-generation sequencing data.

version 1.0

task infer_strandedness {
    input {
        File bam
        File bai
        File gtf
        Int max_retries = 1
        Int memory_gb = 5
    }

    String out_file = basename(bam, ".bam") + ".strandedness.txt"
    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil(((bam_size) * 2) + 10)
 
    command {
        set -euo pipefail

        annotation=~{gtf}
        if [ $(file ~{gtf} | grep -c "gzip compressed data") -eq 1 ]
        then
            gzip -dc ~{gtf} > annotation.gtf
            annotation="annotation.gtf"
        elif [ $(file ~{gtf} | grep -c "bzip2 compressed data") -eq 1 ]
        then
            bzip2 -dc ~{gtf} > annotation.gtf
            annotation="annotation.gtf"
        fi
        
        sort -k1,1 -k4,4n -k5,5n $annotation | bgzip > annotation.gtf.gz
        tabix -p gff annotation.gtf.gz
        mv ~{bai} ~{bam}.bai || true
        ngsderive strandedness ~{bam} -g annotation.gtf.gz > ~{out_file}
        awk 'NR > 1' ~{out_file} | cut -d$'\t' -f5 > strandedness.txt
    }

    runtime {
        disk: disk_size + " GB"
        docker: 'stjudecloud/ngsderive:1.0.2'
        memory: memory_gb + " GB"
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
    }

    String out_file = basename(bam, ".bam") + ".instrument.txt"
    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil(((bam_size) * 1.25) + 10)

    command {
        ngsderive instrument ~{bam} > ~{out_file}
    }

    runtime {
        disk: disk_size + " GB"
        docker: 'stjudecloud/ngsderive:1.0.2'
        maxRetries: max_retries
    }

    output {
        File instrument_file = out_file
    }
}

task read_length {
    input {
        File bam
        File bai
        Int max_retries = 1
        Int memory_gb = 5
    }

    String out_file = basename(bam, ".bam") + ".readlength.txt"
    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil(((bam_size) * 2) + 10)
 
    command {
        set -euo pipefail
        
        mv ~{bai} ~{bam}.bai || true
        ngsderive readlen ~{bam} > ~{out_file}
    }

    runtime {
        disk: disk_size + " GB"
        memory: memory_gb + " GB"
        docker: 'stjudecloud/ngsderive:1.0.2'
        maxRetries: max_retries
    }

    output {
        String computed_read_length = read_string(out_file)
        File read_length_file = out_file
    }

}

task encoding {
    input {
        Array[File] fastqs
        String prefix
        Int num_reads = -1
        Int max_retries = 1
        Int memory_gb = 5
    }

    String out_file = prefix + ".encoding.txt"
    Float fastqs_size = size(fastqs, "GiB")
    Int disk_size = ceil((fastqs_size) + 10)
 
    command {
        ngsderive encoding -n ~{num_reads} ~{sep=' ' fastqs} > ~{out_file}
    }

    runtime {
        disk: disk_size + " GB"
        memory: memory_gb + " GB"
        docker: 'stjudecloud/ngsderive:branch-encoding-1.1.0'
        maxRetries: max_retries
    }

    output {
        File encoding_file = out_file
    }

}
