## # HTSeq
##
## This WDL tool wraps the [htseq](https://github.com/htseq/htseq) tool.
## HTSeq is a Python library for analyzing sequencing data.

version 1.0

task count {
    input {
        File bam
        File gtf
        String outfile_name = basename(bam, ".bam") + ".feature-counts.txt"
        String provided_strandedness = ""
        String inferred_strandedness = ""
        Boolean pos_sorted = true
        Int minaqual = 10
        String feature_type = "exon"
        String idattr = "gene_name"
        String mode = "union"
        Boolean nonunique = false
        Boolean secondary_alignments = false
        Boolean supplementary_alignments = false
        Int added_memory_gb = 0
        Int max_retries = 1
    }

    parameter_meta {
        bam: "Input BAM format file to generate coverage for"
        gtf: "Input genomic features in GTF format to count reads for"
        outfile_name: "Name for the feature count file"
        provided_strandedness: "<Stranded-Reverse/Stranded-Forward/Unstranded> User provided strandedness."
        inferred_strandedness: "<Stranded-Reverse/Stranded-Forward/Unstranded> `ngsderive` derived strandedness."
        pos_sorted: "Is the BAM position sorted (true) or name sorted (false)?"
        minaqual: "Skip all reads with alignment quality lower than the given minimum value"
        feature_type: "Feature type (3rd column in GTF file) to be used, all features of other type are ignored"
        idattr: "GFF attribute to be used as feature ID"
        mode: "<union/intersection-strict/intersection-nonempty> Mode to handle reads overlapping more than one feature"
        nonunique: "Score reads that align to or are assigned to more than one feature?"
        secondary_alignments: "Score secondary alignments (SAM flag 0x100)?"
        supplementary_alignments: "Score supplementary/chimeric alignments (SAM flag 0x800)?"
        added_memory_gb: "Additional RAM to allocate for task. Default RAM is allocated dynamically based on the BAM size."
        max_retries: "Number of times to retry in case of failure"
    }

    String stranded = if (provided_strandedness != "") then 
                        if (provided_strandedness == "Stranded-Reverse") then "reverse" else
                        if (provided_strandedness == "Stranded-Forward") then "yes" else
                        if (provided_strandedness == "Unstranded") then "no"
                         else "unknown-strand" # this will intentionally cause htseq to error.
                                               # You will need to manually specify in this case
                      else 
                        if (inferred_strandedness == "Stranded-Reverse") then "reverse" else
                        if (inferred_strandedness == "Stranded-Forward") then "yes" else 
                        if (inferred_strandedness == "Unstranded") then "no" else
                        if (inferred_strandedness == "Inconclusive") then "no"
                        else "unknown-strand" # this will intentionally cause htseq to error.
                                              # You will need to manually specify in this case

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
            --nonunique ~{if nonunique then "all" else "none"} \
            --secondary-alignments ~{if secondary_alignments then "score" else "ignore"} \
            --supplementary-alignments ~{if supplementary_alignments then "score" else "ignore"} \
            --max-reads-in-buffer 9223372036854776000 \
            ~{bam} \
            ~{gtf} \
            > ~{outfile_name}
    }

    runtime {
        memory: mem_size + " GB"
        disk: disk_size + " GB"
        docker: 'ghcr.io/stjudecloud/htseq:1.0.2'
        maxRetries: max_retries
    }
   
    output {
        File feature_counts = "~{outfile_name}"
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "This WDL tool performs read counting for a set of features in the input BAM file."
    }
}
