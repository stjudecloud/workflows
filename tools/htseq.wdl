## Description:
##
## This WDL tool wraps the htseq tool (https://github.com/simon-anders/htseq).
## HTSeq is a Python library for analyzing sequencing data. 

task count {
    File bam
    File gtf
    String strand = "reverse"
    String outfile = basename(bam, ".bam") + ".counts.txt"

    Int bam_size = size(bam, "GiB")
    Int gtf_size = size(gtf, "GiB")
    Int disk_size = ceil(((bam_size + gtf_size) * 2) + 10)
 
    command {
        htseq-count -f bam -r pos -s ${strand} -m union -i gene_name --secondary-alignments ignore --supplementary-alignments ignore ${bam} ${gtf} > ${outfile}
    }

    runtime {
        memory: "8 GB"
        disk: disk_size + " GB"
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }
   
    output {
       File out = "${outfile}"
    }
}
