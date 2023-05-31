## # HTSeq
##
## This WDL file wraps the [htseq](https://github.com/htseq/htseq) tool.
## HTSeq is a Python library for analyzing sequencing data.

version 1.0

task count {
    meta {
        description: "This WDL task performs read counting for a set of features in the input BAM file."
    }

    parameter_meta {
        bam: "Input BAM format file to generate coverage for"
        gtf: "Input genomic features in gzipped GTF format to count reads for"
        strandedness: {
            description: "Strandedness protocol of the RNA-Seq experiment. Read the [HTSeq-count docs](https://htseq.readthedocs.io/en/release_0.11.1/count.html#cmdoption-htseq-count-s) for an explanation of what this parameter affects."
            choices: [
                "yes",
                "reverse",
                "no"
            ]
        }
        outfile_name: "Name of the feature count file"
        feature_type: "Feature type (3rd column in GTF file) to be used, all features of other type are ignored"
        idattr: "GFF attribute to be used as feature ID"
        mode: {
            description: "Mode to handle reads overlapping more than one feature. `union` is recommended for most use-cases. Read the [HTSeq-count docs](https://htseq.readthedocs.io/en/release_0.11.1/count.html#counting-reads-in-features-with-htseq-count) for an explanation of what this parameter affects."
            choices: [
                "union",
                "intersection-strict",
                "intersection-nonempty"
            ]
        }
        pos_sorted: "Is the BAM position sorted (true) or name sorted (false)?"
        nonunique: "Score reads that align to or are assigned to more than one feature?"
        secondary_alignments: "Score secondary alignments (SAM flag 0x100)?"
        supplementary_alignments: "Score supplementary/chimeric alignments (SAM flag 0x800)?"
        minaqual: "Skip all reads with alignment quality lower than the given minimum value"
        modify_memory_gb: "Add to or subtract from dynamic memory allocation. Default memory is determined by the size of the inputs. Specified in GB."
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
        max_retries: "Number of times to retry in case of failure"
    }

    input {
        File bam
        File gtf
        String strandedness
        String outfile_name = basename(bam, ".bam") + ".feature-counts.txt"
        String feature_type = "exon"
        String idattr = "gene_name"
        String mode = "union"
        Boolean pos_sorted = true
        Boolean nonunique = false
        Boolean secondary_alignments = false
        Boolean supplementary_alignments = false
        Int minaqual = 10
        Int modify_memory_gb = 0
        Int modify_disk_size_gb = 0
        Int max_retries = 1
    }

    Int bam_size = size(bam, "GiB")
    Int gtf_size = size(gtf, "GiB")

    Int memory_gb = bam_size + 5 + modify_memory_gb  # TODO what's the memory usage like if the BAM is name sorted?

    Int disk_size_gb = ceil((bam_size + gtf_size) * 4) + 10 + modify_disk_size_gb  # TODO why would htseq need this much disk?
 
    command <<<
        # 9223372036854776000 == max 64 bit Float
        htseq-count -f bam \
            --max-reads-in-buffer 9223372036854776000 \
            -r ~{if pos_sorted then "pos" else "name"} \
            -s ~{strandedness} \
            -a ~{minaqual} \
            -t ~{feature_type} \
            -m ~{mode} \
            -i ~{idattr} \
            --nonunique ~{if nonunique then "all" else "none"} \
            --secondary-alignments ~{if secondary_alignments then "score" else "ignore"} \
            --supplementary-alignments ~{if supplementary_alignments then "score" else "ignore"} \
            ~{bam} \
            ~{gtf} \
            > ~{outfile_name}
    >>>
   
    output {
        File feature_counts = "~{outfile_name}"
    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        docker: 'quay.io/biocontainers/htseq:2.0.2--py310ha14a713_0'
        maxRetries: max_retries
    }
}
