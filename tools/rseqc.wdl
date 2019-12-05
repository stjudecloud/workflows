## Description:
##
## This WDL tool wraps the RSeQC tool (http://rseqc.sourceforge.net).
## RSeQC is a package for quality control of RNA-seq data.

version 1.0

task infer_experiment {
    input {
        File bam
        Int? sample_size
        Int? map_qual
        File refgene_bed
    }
    Float bam_size = size(bam, "GiB")
    Float refgene_bed_size = size(refgene_bed, "GiB")
    Int disk_size = ceil(((bam_size + refgene_bed_size) * 2) + 10)
 
    command {
        infer_experiment.py -i ${bam} -r ${refgene_bed} ${"-s" + sample_size} ${"-q" + map_qual} > stdout.txt
    }

    runtime {
        disk: disk_size + " GB"
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }

    output {
       String out = read_string("stdout.txt")
    }
}
