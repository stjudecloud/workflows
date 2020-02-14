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
        Int memory_gb = 5
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil(((bam_size) * 2) + 10)
 
    command {
        sort -k1,1 -k4,4n -k5,5n ~{gtf} | bgzip > annotation.gtf.gz
        tabix -p gff annotation.gtf.gz
        mv ~{bai} ~{bam}.bai 
        ngsderive strandedness ~{bam} -g annotation.gtf.gz | awk 'NR > 1' | cut -d$'\t' -f5
    }

    runtime {
        disk: disk_size + " GB"
        memory: memory_gb + " GB"
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
        maxRetries: max_retries
    }

    output {
        String strandedness = read_string(stdout())
        File strandedness_file = stdout()
    }
}
