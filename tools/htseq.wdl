## Description:
##
## This WDL tool wraps the htseq tool (https://github.com/simon-anders/htseq).
## HTSeq is a Python library for analyzing sequencing data.

version 1.0

task count {
    input {
        File bam
        File gtf
        String strand
        String inferred
        String outfile = basename(bam, ".bam") + ".counts.txt"
    }

    String stranded = if (strand != "") then strand else
                 if (inferred == "Stranded-Reverse") then "reverse" else
                 if (inferred == "Stranded-Forward") then "yes" else 
                 if (inferred == "Unstranded") then "no" else
                 "unknown" # this will intentionally cause htseq to error. You will need to manually specify
                           # in this case.

    Float bam_size = size(bam, "GiB")
    Float gtf_size = size(gtf, "GiB")
    Int disk_size = ceil(((bam_size + gtf_size) * 2) + 10)
 
    command {
        htseq-count -f bam -r pos -s ${stranded} -m union -i gene_name --secondary-alignments ignore --supplementary-alignments ignore ${bam} ${gtf} > ${outfile}
    }

    runtime {
        memory: "8 GB"
        disk: disk_size + " GB"
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }
   
    output {
        File out = "${outfile}"
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "This WDL tool performs read counting for a set of features in the input BAM file."
    }

    parameter_meta {
        bam: "Input BAM format file to generate coverage for"
        gtf: "Input genomic features in GTF format to count reads for"
    }
}
