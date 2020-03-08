## Description:
##
## This WDL tool wraps the QualiMap tool (http://qualimap.bioinfo.cipf.es/).
## QualiMap computes metrics to facilitate evaluation of sequencing data. 

version 1.0

task bamqc {
    input {
        File bam
        Int ncpu = 1
        String prefix = basename(bam, ".bam") + "_qualimap_bamqc_results"
        Int max_retries = 1
        Int memory_gb = 8
        Int? disk_size_gb
    }

    Int java_heap_size = ceil(memory_gb * 0.9)
    String out_directory = select_first([prefix, "qualimap_bamqc_results"])
    String out_tar_gz_file = out_directory + ".tar.gz"
    Float bam_size = size(bam, "GiB")
    Int disk_size = select_first([disk_size_gb, ceil((bam_size * 2) + 10)])

    command {
        qualimap bamqc -bam ~{bam} \
            -outdir ~{out_directory} \
            -nt ~{ncpu} \
            --java-mem-size=~{java_heap_size}g
        tar -czf ~{out_tar_gz_file} ~{out_directory}
    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        docker: 'stjudecloud/qualimap:1.0.0-alpha'
        maxRetries: max_retries
    }

    output {
        File results = out_tar_gz_file
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "This WDL tool runs QualiMap's bamqc tool on the input BAM file."
    }

    parameter_meta {
        bam: "Input BAM format file to generate coverage for"
    }
}

task rnaseq {
    input {
        File bam
        File gencode_gtf
        String prefix = basename(bam, ".bam") + "_qualimap_rnaseq_results"
        Int memory_gb = 16
        Int? disk_size_gb
        Int max_retries = 1
        String strand = ""
        String inferred = ""
    }

    String out_directory = select_first([prefix, "qualimap_rnaseq_results"])
    String out_tar_gz_file = out_directory + ".tar.gz"
    String stranded = if (strand != "") then 
                        if (strand == "reverse") then "strand-specific-reverse" else
                        if (strand == "yes") then "strand-specific-forward" else
                        if (strand == "no") then "non-strand-specific"
                        else "unknown-strand"
                      else 
                        if (inferred == "Stranded-Reverse") then "strand-specific-reverse" else
                        if (inferred == "Stranded-Forward") then "strand-specific-forward" else 
                        if (inferred == "Unstranded") then "non-strand-specific" else
                        "unknown-strand" # this will intentionally cause htseq to error. You will need to manually specify
                                         # in this case

    Int java_heap_size = ceil(memory_gb * 0.9)
    Float bam_size = size(bam, "GiB")
    Float gencode_gtf_size = size(gencode_gtf, "GiB")
    Int disk_size = select_first([disk_size_gb, ceil(((bam_size + gencode_gtf_size) * 12) + 10)])
 
    command {
        qualimap rnaseq -bam ~{bam} \
                        -gtf ~{gencode_gtf} \
                        -outdir ~{out_directory} \
                        -oc qualimap_counts.txt \
                        -p ~{stranded} \
                        -pe \
                        --java-mem-size=~{java_heap_size}G
        tar -czf ~{out_tar_gz_file} ~{out_directory}
    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        docker: 'stjudecloud/qualimap:1.0.0-alpha'
        maxRetries: max_retries
    }

    output {
        File results = out_tar_gz_file
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "This WDL tool generates runs QualiMap's rnaseq tool on the input BAM file."
    }

    parameter_meta {
        bam: "Input BAM format file to generate coverage for"
        gencode_gtf: "A GTF format features file containing Gencode features"
        strand: "Strand information for RNA-seq experiments. Options: [strand-specific-forward, strand-specific-reverse, non-strand-specific]"
    }
}
