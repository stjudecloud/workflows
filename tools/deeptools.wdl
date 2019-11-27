## Description: 
##
## This WDL tool wraps the DeepTools tool (https://deeptools.readthedocs.io/en/develop/index.html).
## DeepTools is a suite of Python tools for analysis of high throughput sequencing analysis.

task bamCoverage {
    File bam
    File bai
    String prefix = basename(bam, ".bam")   
 
    command {
        if [ ! -e ${bam}.bai ] 
        then 
           ln -s ${bai} ${bam}.bai
        fi
 
        bamCoverage --bam ${bam} --outFileName ${prefix}.bw --outFileFormat bigwig --numberOfProcessors "max"
    }

    runtime {
        disk: "80 GB"
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }

    output {
        File bigwig = "${prefix}.bw"
    }
} 
