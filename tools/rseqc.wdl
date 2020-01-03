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
    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
        author: "Andrew Frantz"
        email: "andrew.frantz@stjude.org"
        description: "This WDL tool RSeQC's infer_experiment on the input BAM file to infer strandedness information of the underlying RNA-seq experiment."
    }
    parameter_meta {
        bam: "Input BAM format file to generate coverage for"
        refgene_bed: "RefGene features in BED format"
    }
}
