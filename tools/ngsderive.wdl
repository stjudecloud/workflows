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
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil(((bam_size) * 2) + 10)
 
    command {
        sort -k1,1 -k4,4n -k5,5n ~{gtf} | bgzip > annotation.gtf.gz
        tabix -p gff annotation.gtf.gz
        mv ~{bai} ~{bam}.bai 
        ls
        ngsderive strandedness ~{bam} -g annotation.gtf.gz | tail +1 | cut -d$'\t' -f5
    }

    runtime {
        disk: disk_size + " GB"
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }

    output {
        String strandedness = read_string(stdout())
    }
}
