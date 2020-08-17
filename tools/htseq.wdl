## # HTSeq
##
## This WDL tool wraps the [htseq](https://github.com/simon-anders/htseq) tool.
## HTSeq is a Python library for analyzing sequencing data.

version 1.0

task count {
    input {
        File bam
        File gtf
        String provided_strand
        String inferred_strand
        String outfile = basename(bam, ".bam") + ".counts.txt"
        Int added_memory_gb = 20
        Int max_retries = 1
    }

    String stranded = if (provided_strand != "") then 
                        if (provided_strand == "Stranded-Reverse") then "reverse" else
                        if (provided_strand == "Stranded-Forward") then "yes" else
                        if (provided_strand == "Unstranded") then "no"
                        else "unknown-strand"
                      else 
                        if (inferred_strand == "Stranded-Reverse") then "reverse" else
                        if (inferred_strand == "Stranded-Forward") then "yes" else 
                        if (inferred_strand == "Unstranded") then "no" else
                        "unknown-strand" # this will intentionally cause htseq to error. You will need to manually specify
                                         # in this case

    Float bam_size = size(bam, "GiB")
    Float mem_size = bam_size + added_memory_gb
    Float gtf_size = size(gtf, "GiB")
    Int disk_size = ceil(((bam_size + gtf_size) * 4) + 10)
 
    command {
        htseq-count -f bam \
            -r pos \
            -s ~{stranded} \
            -m union \
            -i gene_name \
            --secondary-alignments ignore \
            --supplementary-alignments ignore \
            ~{bam} \
            ~{gtf} \
            > ~{outfile}
    }

    runtime {
        memory: mem_size + " GB"
        disk: disk_size + " GB"
        docker: 'stjudecloud/htseq:1.0.0'
        maxRetries: max_retries
    }
   
    output {
        File out = "~{outfile}"
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "This WDL tool performs read counting for a set of features in the input BAM file."
    }

    parameter_meta {
        bam: "Input BAM format file to generate coverage for"
        gtf: "Input genomic features in GTF format to count reads for"
        added_memory_gb: "Amount of additional memory to add to the bam size"
    }
}
