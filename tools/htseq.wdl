## # HTSeq
##
## This WDL tool wraps the [htseq](https://github.com/htseq/htseq) tool.
## HTSeq is a Python library for analyzing sequencing data.

version 1.0

task count {
    input {
        File bam
        File gtf
        String outfile = basename(bam, ".bam") + ".feature-counts.txt"
        String provided_strandedness = ""
        String inferred_strandedness = ""
        Boolean pos_sorted = true
        Int minaqual = 10
        String feature_type = "exon"
        String idattr = "gene_name"
        String mode = "union"
        String nonunique = "none"
        Boolean secondary_alignments = false
        Boolean supplementary_alignments = false
        Int added_memory_gb = 20
        Int max_retries = 1
    }

    String stranded = if (provided_strandedness != "") then 
                        if (provided_strandedness == "Stranded-Reverse") then "reverse" else
                        if (provided_strandedness == "Stranded-Forward") then "yes" else
                        if (provided_strandedness == "Unstranded") then "no" else
                        "unknown-strand" # this will intentionally cause htseq to error. You will need to manually specify
                                         # in this case
                      else 
                        if (inferred_strandedness == "Stranded-Reverse") then "reverse" else
                        if (inferred_strandedness == "Stranded-Forward") then "yes" else 
                        if (inferred_strandedness == "Unstranded") then "no" else
                        if (inferred_strandedness == "Inconclusive") then "no"
                        else "unknown-strand" # this will intentionally cause htseq to error. You will need to manually specify
                                         # in this case

    Float bam_size = size(bam, "GiB")
    Float mem_size = bam_size + added_memory_gb
    Float gtf_size = size(gtf, "GiB")
    Int disk_size = ceil(((bam_size + gtf_size) * 4) + 10)
 
    command {
        htseq-count -f bam \
            -r ~{if pos_sorted then "pos" else "name"} \
            -s ~{stranded} \
            -a ~{minaqual} \
            -t ~{feature_type} \
            -m ~{mode} \
            -i ~{idattr} \
            --nonunique ~{nonunique} \
            --secondary-alignments ~{if secondary_alignments then "score" else "ignore"} \
            --supplementary-alignments ~{if supplementary_alignments then "score" else "ignore"} \
            --max-reads-in-buffer 9223372036854776000 \
            ~{bam} \
            ~{gtf} \
            > ~{outfile}
    }

    runtime {
        memory: mem_size + " GB"
        disk: disk_size + " GB"
        docker: 'ghcr.io/stjudecloud/htseq:1.0.2'
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
