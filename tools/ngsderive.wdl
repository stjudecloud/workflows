## [Homepage](https://github.com/stjudecloud/ngsderive)

version 1.1

task strandedness {
    meta {
        description: "Derives the experimental strandedness protocol used to generate the input RNA-Seq BAM file. Reports evidence supporting final results."
        outputs: {
            strandedness_file: "TSV file containing the `ngsderive strandedness` report",
            strandedness_string: "The derived strandedness, in string format",
        }
    }

    parameter_meta {
        bam: "Input BAM format file to derive strandedness for"
        bam_index: "BAM index file corresponding to the input BAM"
        gene_model: "Gene model as a GFF/GTF file"
        outfile_name: "Name for the strandedness TSV file"
        split_by_rg: {
            description: "Contain one entry in the output TSV per read group, in addition to an `overall` entry",
            group: "Common",
        }
        min_reads_per_gene: {
            description: "Filter any genes that don't have at least `min_reads_per_gene` reads mapping to them",
            group: "Common",
        }
        num_genes: {
            description: "How many genes to sample",
            group: "Common",
        }
        min_mapq: {
            description: "Minimum MAPQ to consider for supporting reads",
            group: "Common",
        }
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        File bam_index
        File gene_model
        String outfile_name = basename(bam, ".bam") + ".strandedness.tsv"
        Boolean split_by_rg = false
        Int min_reads_per_gene = 10
        Int num_genes = 1000
        Int min_mapq = 30
        Int modify_disk_size_gb = 0
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        # localize BAM, BAI, and GFF to CWD
        CWD_BAM=~{basename(bam)}
        ln -s "~{bam}" "$CWD_BAM"
        ln -s "~{bam_index}" "$CWD_BAM".bai
        CWD_GFF=~{basename(gene_model)}
        ln -s "~{gene_model}" "$CWD_GFF"

        ngsderive strandedness --verbose \
            ~{if split_by_rg then "--split-by-rg" else ""} \
            -m ~{min_reads_per_gene} \
            -n ~{num_genes} \
            -q ~{min_mapq} \
            -g "$CWD_GFF" \
            "$CWD_BAM" \
            > "~{outfile_name}"

        if ~{split_by_rg}; then
            echo "N/A" > strandedness.txt
        else
            awk 'NR > 1' "~{outfile_name}" | cut -f 5 > strandedness.txt
        fi

        rm "$CWD_BAM" "$CWD_BAM".bai
    >>>

    output {
        File strandedness_file = outfile_name
        String strandedness_string = read_string("strandedness.txt")
    }

    runtime {
        memory: "4 GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/ngsderive:3.3.2--pyhdfd78af_0"
        maxRetries: 1
    }
}

task instrument {
    meta {
        description: "Derives the instrument used to sequence the input BAM file. Reports evidence supporting final results."
        outputs: {
            instrument_file: "TSV file containing the `ngsderive isntrument` report for the input BAM file",
            instrument_string: "The derived instrument, in string format",
        }
    }

    parameter_meta {
        bam: "Input BAM format file to derive instrument for"
        outfile_name: "Name for the instrument TSV file"
        num_reads: {
            description: "How many reads to analyze from the start of the file. Any n < 1 to parse whole file.",
            group: "Common",
        }
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        String outfile_name = basename(bam, ".bam") + ".instrument.tsv"
        Int num_reads = 10000
        Int modify_disk_size_gb = 0
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        ngsderive instrument --verbose \
            -n ~{num_reads} \
            "~{bam}" \
            > "~{outfile_name}"

        awk 'NR > 1' "~{outfile_name}" | cut -f 2 > instrument.txt
    >>>

    output {
        File instrument_file = outfile_name
        String instrument_string = read_string("instrument.txt")
    }

    runtime {
        memory: "4 GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/ngsderive:3.3.2--pyhdfd78af_0"
        maxRetries: 1
    }
}

task read_length {
    meta {
        description: "Derives the original experimental read length of the input BAM. Reports evidence supporting final results."
        outputs: {
            read_length_file: "TSV file containing the `ngsderive readlen` report for the input BAM file",
        }
    }

    parameter_meta {
        bam: "Input BAM format file to derive read length for"
        bam_index: "BAM index file corresponding to the input BAM"
        outfile_name: "Name for the readlen TSV file"
        majority_vote_cutoff: {
            description: "To call a majority readlen, the maximum read length must have at least `majority-vote-cutoff`% reads in support",
            group: "Common",
        }
        num_reads: {
            description: "How many reads to analyze from the start of the file. Any n < 1 to parse whole file.",
            group: "Common",
        }
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        File bam_index
        String outfile_name = basename(bam, ".bam") + ".readlength.tsv"
        Float majority_vote_cutoff = 0.7
        Int num_reads = -1
        Int modify_disk_size_gb = 0
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        # localize BAM and BAI to CWD
        CWD_BAM=~{basename(bam)}
        ln -s "~{bam}" "$CWD_BAM"
        ln -s "~{bam_index}" "$CWD_BAM".bai

        ngsderive readlen --verbose \
            -c ~{majority_vote_cutoff} \
            -n ~{num_reads} \
            "$CWD_BAM" \
            > "~{outfile_name}"

        rm "$CWD_BAM" "$CWD_BAM".bai
    >>>

    output {
        File read_length_file = outfile_name
    }

    runtime {
        memory: "4 GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/ngsderive:3.3.2--pyhdfd78af_0"
        maxRetries: 1
    }
}

task encoding {
    meta {
        description: "Derives the encoding of the input NGS file(s). Reports evidence supporting final results."
        outputs: {
            encoding_file: "TSV file containing the `ngsderive encoding` report for all input files",
        }
    }

    parameter_meta {
        ngs_files: "An array of FASTQs and/or BAMs for which to derive encoding"
        outfile_name: "Name for the encoding TSV file"
        num_reads: {
            description: "How many reads to analyze from the start of the file(s). Any n < 1 to parse whole file(s).",
            group: "Common",
        }
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        Array[File] ngs_files
        String outfile_name
        Int num_reads = 1000000
        Int modify_disk_size_gb = 0
    }

    Float files_size = size(ngs_files, "GiB")
    Int disk_size_gb = ceil(files_size) + 10 + modify_disk_size_gb

    command <<<
        ngsderive encoding --verbose \
            -n ~{num_reads} \
            ~{sep(" ", squote(ngs_files))} \
            > "~{outfile_name}"
    >>>

    output {
        File encoding_file = outfile_name
    }

    runtime {
        memory: "4 GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/ngsderive:3.3.2--pyhdfd78af_0"
        maxRetries: 1
    }
}

task junction_annotation {
    meta {
        description: "Annotates junctions found in an RNA-Seq BAM as known, novel, or partially novel"
        external_help: "https://stjudecloud.github.io/ngsderive/subcommands/junction_annotation/"
        outputs: {
            junction_summary: "TSV file containing the `ngsderive junction-annotation` summary",
            junctions: "TSV file containing a detailed list of annotated junctions",
        }
    }

    parameter_meta {
        bam: "Input BAM format file to annotate junctions for"
        bam_index: "BAM index file corresponding to the input BAM"
        gene_model: "Gene model as a GFF/GTF file"
        prefix: "Prefix for the summary TSV and junction files. The extensions `.junction_summary.tsv` and `.junctions.tsv` will be added."
        min_intron: {
            description: "Minimum size of intron to be considered a splice",
            group: "Common",
        }
        min_mapq: {
            description: "Minimum MAPQ to consider for supporting reads",
            group: "Common",
        }
        min_reads: {
            description: "Filter any junctions that don't have at least `min_reads` reads supporting them",
            group: "Common",
        }
        fuzzy_junction_match_range: {
            description: "Consider found splices within `+-k` bases of a known splice event annotated",
            group: "Common",
        }
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        File bam_index
        File gene_model
        String prefix = basename(bam, ".bam")
        Int min_intron = 50
        Int min_mapq = 30
        Int min_reads = 2
        Int fuzzy_junction_match_range = 0
        Int modify_disk_size_gb = 0
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        # localize BAM and BAI to CWD
        CWD_BAM=~{basename(bam)}
        ln -s "~{bam}" "$CWD_BAM"
        ln -s "~{bam_index}" "$CWD_BAM".bai

        ngsderive junction-annotation --verbose \
            -g "~{gene_model}" \
            -i ~{min_intron} \
            -q ~{min_mapq} \
            -m ~{min_reads} \
            -k ~{fuzzy_junction_match_range} \
            -o "~{prefix}.junction_summary.tsv" \
            "$CWD_BAM"

        # junction-annotation accepts multiple BAMs, and allows for
        # renaming the cohort level summary report, but not the BAM
        # level junctions file. So we rename it here to match with prefix.
        mv "$(basename "~{bam}.junctions.tsv")" "~{prefix}.junctions.tsv"
        gzip "~{prefix}.junctions.tsv"

        rm "$CWD_BAM" "$CWD_BAM".bai
    >>>

    output {
        File junction_summary = "~{prefix}.junction_summary.tsv"
        File junctions = "~{prefix}.junctions.tsv.gz"
    }

    runtime {
        memory: "56 GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/ngsderive:3.3.2--pyhdfd78af_0"
        maxRetries: 1
    }
}

task endedness {
    meta {
        description: "Derives the endedness of the input BAM file. Reports evidence for final result."
        outputs: {
            endedness_file: "TSV file containing the `ngsderive endedness` report",
        }
    }

    parameter_meta {
        bam: "Input BAM format file to derive endedness from"
        outfile_name: "Name for the endedness TSV file"
        lenient: {
            description: "Return a zero exit code on unknown results",
            group: "Common",
        }
        calc_rpt: {
            description: "Calculate and output Reads-Per-Template.",
            warning: "This will produce a more sophisticated estimate for endedness, but uses substantially more memory (can reach up to 200% of BAM size in memory consumption for some inputs).",
            group: "Common",
        }
        round_rpt: {
            description: "Round RPT to the nearest INT before comparing to expected values. Appropriate if using `--num-reads` > 0.",
            group: "Common",
        }
        split_by_rg: {
            description: "Contain one entry per read group",
            group: "Common",
        }
        paired_deviance: {
            description: "Distance from 0.5 split between number of f+l- reads and f-l+ reads allowed to be called 'Paired-End'.",
            warning: "Default of `0.0` only appropriate if the whole file is being processed.",
            group: "Common",
        }
        num_reads: {
            description: "How many reads to analyze from the start of the file. Any n < 1 to parse whole file.",
            group: "Common",
        }
        modify_memory_gb: "Add to or subtract from dynamic memory allocation. Default memory is determined by value of `calc_rpt` and size of the BAM. Specified in GB."
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        String outfile_name = basename(bam, ".bam") + ".endedness.tsv"
        Boolean lenient = false
        Boolean calc_rpt = false
        Boolean round_rpt = false
        Boolean split_by_rg = false
        Float paired_deviance = 0.0
        Int num_reads = -1
        Int modify_memory_gb = 0
        Int modify_disk_size_gb = 0
    }

    Float bam_size = size(bam, "GiB")
    Int memory_gb = (
        if calc_rpt
        then (
            ceil(bam_size * 2.5) + 4 + modify_memory_gb
        )
        else 4
    )
    Int disk_size_gb = ceil(bam_size) + 10 + modify_disk_size_gb

    command <<<
        ngsderive endedness --verbose \
            ~{if lenient then "--lenient" else ""} \
            ~{if calc_rpt then "-r" else ""} \
            ~{if round_rpt then "--round-rpt" else ""} \
            ~{if split_by_rg then "--split-by-rg" else ""} \
            --paired-deviance ~{paired_deviance} \
            -n ~{num_reads} \
            "~{bam}" \
            > "~{outfile_name}"
    >>>

    output {
        File endedness_file = outfile_name
    }

    runtime {
        memory: "~{memory_gb} GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/ngsderive:3.3.2--pyhdfd78af_0"
        maxRetries: 1
    }
}
