## # HTSeq
##
## This WDL file wraps the [htseq](https://github.com/htseq/htseq) tool.
## HTSeq is a Python library for analyzing sequencing data.

version 1.0

task count {
    input {
        File bam
        File gtf
        String provided_strandedness
        String inferred_strandedness = ""
        String outfile = basename(bam, ".bam") + ".feature-counts.txt"
        Int added_memory_gb = 20
        Int max_retries = 1
    }

    String stranded = if (provided_strandedness != "") then 
                        if (provided_strandedness == "Stranded-Reverse") then "reverse" else
                        if (provided_strandedness == "Stranded-Forward") then "yes" else
                        if (provided_strandedness == "Unstranded") then "no"
                        else "unknown-strand"
                      else 
                        if (inferred_strandedness == "Stranded-Reverse") then "reverse" else
                        if (inferred_strandedness == "Stranded-Forward") then "yes" else 
                        if (inferred_strandedness == "Unstranded") then "no" else
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
            --max-reads-in-buffer 9223372036854776000 \
            ~{bam} \
            ~{gtf} \
            > ~{outfile}
    }

    runtime {
        memory: mem_size + " GB"
        disk: disk_size + " GB"
        docker: 'quay.io/biocontainers/htseq:2.0.2--py310ha14a713_0'
        maxRetries: max_retries
    }
   
    output {
        File out = "~{outfile}"
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "This WDL task performs read counting for a set of features in the input BAM file."
    }

    parameter_meta {
        bam: "Input BAM format file to generate coverage for"
        gtf: "Input genomic features in GTF format to count reads for"
        added_memory_gb: "Amount of additional memory to add to the bam size"
    }
}
