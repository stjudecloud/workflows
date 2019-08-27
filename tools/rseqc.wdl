## Description: 
##
## This WDL tool wraps the RSeQC tool (http://rseqc.sourceforge.net).
## RSeQC is a package for quality control of RNA-seq data.

task infer_experiment {
    File bam
    Int? sample_size
    Int? map_qual
    File refgene_bed 
 
    command {
        infer_experiment.py -i ${bam} -r ${refgene_bed} ${"-s" + sample_size} ${"-q" + map_qual} > stdout.txt 
    }

    runtime {
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }

    output {
       String out = read_string("stdout.txt")
    }
}
