## Description: 
##
## This WDL tool wraps the htseq tool (https://github.com/simon-anders/htseq).
## HTSeq is a Python library for analyzing sequencing data. 

task count {
    File bam
    File gtf
    String strand = "reverse"
    String outfile = basename(bam, ".bam") + ".counts.txt"
 
    command {
        htseq-count -f bam -r pos -s ${strand} -m union -i gene_name --secondary-alignments ignore --supplementary-alignments ignore ${bam} ${gtf} > ${outfile}
    }

    runtime {
        memory: "8G"
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }
   
    output {
       File out = "${outfile}"
    }
    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
        author: "Andrew Frantz"
        email: "andrew.frantz@stjude.org"
        description: "This WDL tool performs read counting for a set of features in the input BAM file."
    }
    parameter_meta {
        bam: "Input BAM format file to generate coverage for"
        gtf: "Input genomic features in GTF format to count reads for"
    }
}
