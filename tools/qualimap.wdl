## # QualiMap
##
## This WDL file wraps the [QualiMap](http://qualimap.bioinfo.cipf.es/) tool.
## QualiMap computes metrics to facilitate evaluation of sequencing data. 

version 1.0

task bamqc {
    input {
        File bam
        Int ncpu = 1
        Int max_retries = 1
        Int memory_gb = 32
        Int? disk_size_gb
        Boolean use_all_cores = false
    }

    String out_directory = basename(bam, ".bam") + '.qualimap_bamqc_results'
    String out_tar_gz = out_directory + ".tar.gz"

    Int java_heap_size = ceil(memory_gb * 0.9)

    Float bam_size = size(bam, "GiB")
    Int disk_size = select_first([disk_size_gb, ceil((bam_size * 2) + 10)])

    command {
        set -euo pipefail
        
        n_cores=~{ncpu}
        if [ "~{use_all_cores}" = "true" ]; then
            n_cores=$(nproc)
        fi
        
        qualimap bamqc -bam ~{bam} \
            -outdir ~{out_directory} \
            -nt "$n_cores" \
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
        docker: 'quay.io/biocontainers/qualimap:2.2.2d--hdfd78af_2'
        maxRetries: max_retries
    }

    output {
        File results = out_tar_gz
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "[Deprecated] This WDL task runs QualiMap's bamqc tool on the input BAM file. It has been deprecated in our pipeline due to memory leak issues. Use at your own risk, for some samples can consume over 1TB of RAM."
    }

    parameter_meta {
        bam: "Input BAM format file to run qualimap bamqc on"
    }
}

task rnaseq {
    input {
        File bam
        File gtf
        String prefix = basename(bam, ".bam")
        Boolean name_sorted = false
        Boolean paired_end = true
        Int memory_gb = 16
        Int? disk_size_gb
        Int max_retries = 1
    }

    parameter_meta {
        bam: "Input BAM format file to run qualimap rnaseq on"
        gtf: "GTF features file"
        prefix: "Prefix for the results directory and output tarball. The extension `.qualimap_rnaseq_results.tar.gz` will be added."
        name_sorted: "Is the BAM name sorted?"
        paired_end: "Is the BAM paired end?"
        memory_gb: "RAM to allocate for task"
        disk_size_gb: "Disk space to allocate for task. Default is determined dynamically based on BAM and GTF sizes."
        max_retries: "Number of times to retry in case of failure"
    }

    String out_directory = prefix + ".qualimap_rnaseq_results"
    String out_tar_gz = out_directory + ".tar.gz"
    String name_sorted_arg = if (name_sorted) then "-s" else ""
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
                        ~{name_sorted_arg} \
                        ~{paired_end_arg} \
                        --java-mem-size=~{java_heap_size}G
        rm "$gtf_name"
        
        # Check if qualimap succeeded
        if [ ! -d "~{out_directory}/raw_data_qualimapReport/" ]; then
            exit 42
        fi
        
        tar -czf ~{out_tar_gz} ~{out_directory}
    >>>

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        docker: 'quay.io/biocontainers/qualimap:2.2.2d--hdfd78af_2'
        maxRetries: max_retries
    }

    output {
        File raw_summary = "~{out_directory}/rnaseq_qc_results.txt"
        File raw_coverage = "~{out_directory}/raw_data_qualimapReport/coverage_profile_along_genes_(total).txt"
        File results = out_tar_gz
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "This WDL task generates runs QualiMap's rnaseq tool on the input BAM file. Note that we don't expose the `-p` parameter. This is used to set strandedness protocol of the sample, however in practice it only disables certain calculations. We do not expose the parameter so that the full suite of calculations is always performed."
    }
}
