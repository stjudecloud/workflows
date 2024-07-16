## [Homepage](https://github.com/htseq/htseq)

version 1.1

task count {
    meta {
        description: "Performs read counting for a set of features in the input BAM file"
        outputs: {
            feature_counts: "A two column headerless TSV file. First column is feature names and second column is counts."
        }
    }

    parameter_meta {
        bam: "Input BAM format file to generate coverage for"
        gtf: "Input genomic features in gzipped GTF format to count reads for"
        strandedness: {
            description: "Strandedness protocol of the RNA-Seq experiment",
            external_help: "https://htseq.readthedocs.io/en/latest/htseqcount.html#cmdoption-htseq-count-s",
            choices: [
                "yes",
                "reverse",
                "no"
            ]
        }
        prefix: "Prefix for the feature counts file. The extension `.feature-counts.txt` will be added."
        feature_type: {
            description: "Feature type (3rd column in GTF file) to be used, all features of other type are ignored",
            common: true
        }
        idattr: {
            description: "GFF attribute to be used as feature ID",
            common: true
        }
        mode: {
            description: "Mode to handle reads overlapping more than one feature. `union` is recommended for most use-cases.",
            external_help: "https://htseq.readthedocs.io/en/latest/htseqcount.html#htseq-count-counting-reads-within-features",
            choices: [
                "union",
                "intersection-strict",
                "intersection-nonempty"
            ]
        }
        include_custom_header: {
            description: "Include a custom header for the output file? This is not an official feature of HTSeq. If true, the first line of the output file will be `feature\t~{prefix}`. This may break downstream tools that expect the typical headerless HTSeq output format.",
            common: true
        }
        pos_sorted: {
            description: "Is the BAM position sorted (true) or name sorted (false)?",
            common: true
        }
        nonunique: {
            description: "Score reads that align to or are assigned to more than one feature?",
            common: true
        }
        secondary_alignments: {
            description: "Score secondary alignments (SAM flag 0x100)?",
            common: true
        }
        supplementary_alignments: {
            description: "Score supplementary/chimeric alignments (SAM flag 0x800)?",
            common: true
        }
        minaqual: {
            description: "Skip all reads with alignment quality lower than the given minimum value",
            common: true
        }
        modify_memory_gb: "Add to or subtract from dynamic memory allocation. Default memory is determined by the size of the inputs. Specified in GB."
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        File gtf
        String strandedness
        String prefix = basename(bam, ".bam")
        String feature_type = "exon"
        String idattr = "gene_name"
        String mode = "union"
        Boolean include_custom_header = false
        Boolean pos_sorted = true
        Boolean nonunique = false
        Boolean secondary_alignments = false
        Boolean supplementary_alignments = false
        Int minaqual = 10
        Int modify_memory_gb = 0
        Int modify_disk_size_gb = 0
    }

    String outfile_name = prefix + ".feature-counts.txt"

    Float bam_size = size(bam, "GiB")
    Float gtf_size = size(gtf, "GiB")

    Int memory_gb = ceil(bam_size) + 5 + modify_memory_gb

    Int disk_size_gb = ceil((bam_size + gtf_size) * 4) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        if ~{include_custom_header}; then
            echo -e "feature\t~{prefix}" > ~{outfile_name}
        else
            true > ~{outfile_name}  # ensure file is empty
        fi

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
            --supplementary-alignments ~{(
                if supplementary_alignments
                then "score"
                else "ignore"
            )} \
            ~{bam} \
            ~{gtf} \
            >> ~{outfile_name}
    >>>

    output {
        File feature_counts = "~{outfile_name}"
    }

    runtime {
        memory: "~{memory_gb} GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/htseq:2.0.5--py310h5aa3a86_0"
        maxRetries: 1
    }
}

task calc_tpm {
    meta {
        description: "Given a gene counts file and a gene lengths file, calculate Transcripts Per Million (TPM)"
        outputs: {
            tpm_file: "Transcripts Per Million (TPM) file. A two column headered TSV file."
        }
    }

    parameter_meta {
        counts: "A two column headerless TSV file with gene names in the first column and counts (as integers) in the second column. Entries starting with '__' will be discarded. Can be generated with the `count` task."
        gene_lengths: "A two column headered TSV file with gene names (matching those in the `counts` file) in the first column and feature lengths (as integers) in the second column. Can be generated with the `calc_gene_lengths` task in `util.wdl`."
        prefix: "Prefix for the TPM file. The extension `.TPM.txt` will be added."
    }

    input {
        File counts
        File gene_lengths
        String prefix = basename(counts, ".feature-counts.txt")
    }

    String outfile_name = prefix + ".TPM.txt"

    #@ except: LineWidth
    command <<<
        COUNTS="~{counts}" GENE_LENGTHS="~{gene_lengths}" OUTFILE="~{outfile_name}" python3 - <<END
        import os  # lint-check: ignore

        counts_file = open(os.environ['COUNTS'], 'r')
        counts = {}
        for line in counts_file:
            gene, count = line.split('\t')
            if gene[0:2] == '__':
                break
            counts[gene.strip()] = int(count.strip())
        counts_file.close()

        lengths_file = open(os.environ['GENE_LENGTHS'], 'r')
        rpks = {}  # Reads Per Kilobase
        tot_rpk = 0
        lengths_file.readline()  # discard header
        for line in lengths_file:
            gene, length = line.split('\t')
            rpk = counts[gene.strip()] / int(length.strip()) * 1000
            tot_rpk += rpk
            rpks[gene.strip()] = rpk
        lengths_file.close()

        scaling_factor = tot_rpk / 1000000

        sample_name = '.'.join(os.environ['OUTFILE'].split('.')[:-2])  # equivalent to ~{prefix}
        outfile = open(os.environ['OUTFILE'], 'w')
        print(f"feature\t{sample_name}", file=outfile)
        for gene, rpk in sorted(rpks.items()):
            tpm = rpk / scaling_factor
            print(f"{gene}\t{tpm:.3f}", file=outfile)
        outfile.close()
        END
    >>>

    output {
        File tpm_file = "~{outfile_name}"
    }

    runtime {
        memory: "4 GB"
        disks: "10 GB"
        container: "ghcr.io/stjudecloud/util:1.3.0"
        maxRetries: 1
    }
}
