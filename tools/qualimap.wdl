## [Homepage](http://qualimap.bioinfo.cipf.es/)
#
# SPDX-License-Identifier: MIT
# Copyright St. Jude Children's Research Hospital
version 1.1

task rnaseq {
    meta {
        description: "Generates runs QualiMap's rnaseq tool on the input BAM file. Note that we don't expose the `-p` parameter. This is used to set strandedness protocol of the sample, however in practice it only disables certain calculations. We do not expose the parameter so that the full suite of calculations is always performed."
        outputs: {
            raw_summary: "Raw text summary of QualiMap's results. Can be parsed by MultiQC."
            raw_coverage: "Raw text of QualiMap's coverage analysis results. Can be parsed by MultiQC."
            results: "Gzipped tar archive of all QualiMap output files"
        }
    }

    parameter_meta {
        bam: "Input BAM format file to run qualimap rnaseq on"
        gtf: "GTF features file. Gzipped or uncompressed."
        prefix: "Prefix for the results directory and output tarball. The extension `.qualimap_rnaseq_results.tar.gz` will be added."
        memory_gb: "RAM to allocate for task"
        name_sorted: {
            description: "Is the BAM name sorted? QualiMap has an inefficient sorting algorithm. In order to save resources we recommend collating your input BAM before QualiMap and setting this parameter to true."
            common: true
        }
        paired_end: {
            description: "Is the BAM paired end?"
            common: true
        }
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        File gtf
        String prefix = basename(bam, ".bam") + ".qualimap_rnaseq_results"
        Boolean name_sorted = false
        Boolean paired_end = true
        Int memory_gb = 16
        Int modify_disk_size_gb = 0
    }

    String out_tar_gz = prefix + ".tar.gz"
    String name_sorted_arg = if (name_sorted) then "-s" else ""
    String paired_end_arg = if (paired_end) then "-pe" else ""

    Int java_heap_size = ceil(memory_gb * 0.9)
    Float bam_size = size(bam, "GiB")
    Float gtf_size = size(gtf, "GiB")

    # Qualimap has an inefficient name sorting algorithm and will
    # use an excessive amount of storage.
    Int disk_size_gb = (
        (
            if name_sorted
            then ceil(bam_size + gtf_size + 15)
            else ceil(((bam_size + gtf_size) * 12) + 10)
        ) + modify_disk_size_gb
    )

    command <<<
        set -euo pipefail

        orig=~{gtf}
        gtf_name=$(basename "${orig%.gz}")
        gunzip -c ~{gtf} > "$gtf_name" || ln -sf ~{gtf} "$gtf_name"

        # '-oc qualimap_counts.txt' puts the file in '-outdir'
        qualimap rnaseq -bam ~{bam} \
                        -oc qualimap_counts.txt \
                        -gtf "$gtf_name" \
                        -outdir ~{prefix} \
                        ~{name_sorted_arg} \
                        ~{paired_end_arg} \
                        --java-mem-size=~{java_heap_size}G
        rm "$gtf_name"

        tar -czf ~{out_tar_gz} ~{prefix}
    >>>

    output {
        File raw_summary = "~{prefix}/rnaseq_qc_results.txt"
        File raw_coverage = "~{prefix}/raw_data_qualimapReport/coverage_profile_along_genes_(total).txt"
        File results = out_tar_gz
    }

    runtime {
        memory: "~{memory_gb} GB"
        disk: "~{disk_size_gb} GB"
        container: 'quay.io/biocontainers/qualimap:2.2.2d--hdfd78af_2'
        maxRetries: 1
    }
}

task bamqc {
    meta {
        description: "**[Deprecated]** This WDL task runs QualiMap's bamqc tool on the input BAM file. This task has been deprecated due to memory leak issues. Use at your own risk, for some samples can consume over 1TB of RAM."
        deprecated: true
    }

    parameter_meta {
        bam: "Input BAM format file to run qualimap bamqc on"
    }

    input {
        File bam
        String prefix = basename(bam, ".bam")
        Boolean use_all_cores = false
        Int ncpu = 1
        Int memory_gb = 32
        Int modify_disk_size_gb = 0
        Int max_retries = 1
    }

    String out_directory = prefix + '.qualimap_bamqc_results'
    String out_tar_gz = out_directory + ".tar.gz"

    Int java_heap_size = ceil(memory_gb * 0.9)

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi

        qualimap bamqc -bam ~{bam} \
            -outdir ~{out_directory} \
            -nt "$n_cores" \
            -nw 400 \
            --java-mem-size=~{java_heap_size}g

        # Check if qualimap succeeded
        if [ ! -d "~{out_directory}/raw_data_qualimapReport/" ]; then
            exit 42
        fi

        tar -czf ~{out_tar_gz} ~{out_directory}
    >>>

    output {
        File results = out_tar_gz
    }

    runtime {
        cpu: ncpu
        memory: "~{memory_gb} GB"
        disk: "~{disk_size_gb} GB"
        container: 'quay.io/biocontainers/qualimap:2.2.2d--hdfd78af_2'
        maxRetries: max_retries
    }
}
