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
        File gtf
        Int max_retries = 1
    }

    String out_file_name = basename(bam, ".bam") + ".strandedness.txt"
    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil(((bam_size) * 2) + 10)
 
    command {
        sort -k1,1 -k4,4n -k5,5n ~{gtf} | bgzip > annotation.gtf.gz
        tabix -p gff annotation.gtf.gz
        mv ~{bai} ~{bam}.bai 
        ngsderive strandedness ~{bam} -g annotation.gtf.gz \
            | awk 'NR > 1' | cut -d$'\t' -f5 > ${out_file_name}
    }

    File out_file = out_file_name

    runtime {
        disk: disk_size + " GB"
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
        maxRetries: max_retries
    }

    output {
        String strandedness = read_string(out_file)
        File strandedness_file = out_file
    }
}
