## # QualiMap
##
## This WDL tool wraps the [QualiMap](http://qualimap.bioinfo.cipf.es/) tool.
## QualiMap computes metrics to facilitate evaluation of sequencing data. 

version 1.0

task bamqc {
    input {
        File bam
        Int ncpu = 1
        Int max_retries = 1
        Int memory_gb = 32
        Int? disk_size_gb
    }

    String out_directory = basename(bam, ".bam") + '.qualimap_bamqc_results'
    String out_tar_gz = out_directory + ".tar.gz"

    Int java_heap_size = ceil(memory_gb * 0.9)
    Float bam_size = size(bam, "GiB")
    Int disk_size = select_first([disk_size_gb, ceil((bam_size * 2) + 10)])

    command {
        set -euo pipefail
        
        qualimap bamqc -bam ~{bam} \
            -outdir ~{out_directory} \
            -nt ~{ncpu} \
            -nw 400 \
            --java-mem-size=~{java_heap_size}g
        
        # Check if qualimap succeeded
        if [ ! -d "~{out_directory}/raw_data_qualimapReport/" ]; then
            exit 1
        fi

        tar -czf ~{out_tar_gz} ~{out_directory}
    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        cpu: ncpu
        docker: 'ghcr.io/stjudecloud/qualimap:1.0.4'
        maxRetries: max_retries
    }

    output {
        File results = out_tar_gz
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
        File gtf
        Int memory_gb = 16
        Int? disk_size_gb
        Int max_retries = 1
        Boolean paired_end = false
        String provided_strandedness = ""
        String inferred_strandedness = ""
    }

    String out_directory = basename(bam, ".bam") + ".qualimap_rnaseq_results"
    String out_tar_gz = out_directory + ".tar.gz"
    String stranded = if (provided_strandedness != "") then 
                        if (provided_strandedness == "Stranded-Reverse") then "strand-specific-reverse" else
                        if (provided_strandedness == "Stranded-Forward") then "strand-specific-forward" else
                        if (provided_strandedness == "Unstranded") then "non-strand-specific"
                        else "unknown-strand"
                      else 
                        if (inferred_strandedness == "Stranded-Reverse") then "strand-specific-reverse" else
                        if (inferred_strandedness == "Stranded-Forward") then "strand-specific-forward" else 
                        if (inferred_strandedness == "Unstranded") then "non-strand-specific" else
                        "unknown-strand" # this will intentionally cause qualimap to error. You will need to manually specify
                                         # in this case
    String paired_end_arg = if (paired_end) then "-pe" else ""

    Int java_heap_size = ceil(memory_gb * 0.9)
    Float bam_size = size(bam, "GiB")
    Float gtf_size = size(gtf, "GiB")
    Int disk_size = select_first([disk_size_gb, ceil(((bam_size + gtf_size) * 12) + 10)])
 
    command <<<
        set -euo pipefail

        orig=~{gtf}
        gtf_name=$(basename "${orig%.gz}")
        gunzip -c ~{gtf} > "$gtf_name" || cp ~{gtf} "$gtf_name"
        
        qualimap rnaseq -bam ~{bam} \
                        -gtf "$gtf_name" \
                        -outdir ~{out_directory} \
                        -oc qualimap_counts.txt \
                        -p ~{stranded} \
                        ~{paired_end_arg} \
                        --java-mem-size=~{java_heap_size}G
        rm "$gtf_name"
        
        # Check if qualimap succeeded
        if [ ! -d "~{out_directory}/raw_data_qualimapReport/" ]; then
            exit 1
        fi
        
        tar -czf ~{out_tar_gz} ~{out_directory}
    >>>

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        docker: 'ghcr.io/stjudecloud/qualimap:1.0.4'
        maxRetries: max_retries
    }

    output {
        File results = out_tar_gz
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "This WDL tool generates runs QualiMap's rnaseq tool on the input BAM file."
    }

    parameter_meta {
        bam: "Input BAM format file to generate coverage for"
        gtf: "A GTF format features file"
        provided_strandedness: "Strand information for RNA-seq experiments. Options: [Stranded-Reverse, Stranded-Forward, Unstranded]"
    }
}
