## # ngsderive
##
## This WDL file wraps the [ngsderive package](https://github.com/stjudecloud/ngsderive)

version 1.1

task strandedness {
    input {
        File bam
        File bam_index
        File gtf
        String outfile_name = basename(bam, ".bam") + ".strandedness.tsv"
        Boolean split_by_rg = false
        Int min_reads_per_gene = 10
        Int num_genes = 1000
        Int min_mapq = 30
        Int memory_gb = 4
        Int modify_disk_size_gb = 0
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size) + 10 + modify_disk_size_gb
 
    command <<<
        set -euo pipefail

        # localize BAM and BAI to CWD
        # some backends prevent writing to the inputs directories
        # to accomodate this, create symlinks in CWD
        CWD_BAM=~{basename(bam)}
        ln -s ~{bam} "$CWD_BAM"
        ln -s ~{bam_index} "$CWD_BAM".bai

        ngsderive strandedness --verbose \
            ~{if split_by_rg then "--split-by-rg" else ""} \
            -m ~{min_reads_per_gene} \
            -n ~{num_genes} \
            -q ~{min_mapq} \
            -g ~{gtf} \
            "$CWD_BAM" \
            > ~{outfile_name}

        awk 'NR > 1' ~{outfile_name} | cut -d$'\t' -f5 > strandedness.txt

        rm "$CWD_BAM" "$CWD_BAM".bai
    >>>

    output {
        String strandedness = read_string("strandedness.txt")
        File strandedness_file = outfile_name
    }

    runtime {
        memory: "~{memory_gb} GB"
        disk: "~{disk_size_gb} GB"
        docker: 'quay.io/biocontainers/ngsderive:3.3.2--pyhdfd78af_0'
        maxRetries: max_retries
    }
}

task instrument {
    input {
        File bam
        String outfile_name = basename(bam, ".bam") + ".instrument.tsv"
        Int num_reads = 10000
        Int memory_gb = 4
        Int modify_disk_size_gb = 0
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size) + 10 + modify_disk_size_gb

    command <<<
        ngsderive instrument --verbose \
            -n ~{num_reads} \
            ~{bam} \
            > ~{outfile_name}
    >>>

    output {
        File instrument_file = outfile_name
    }

    runtime {
        memory: "~{memory_gb} GB"
        disk: "~{disk_size_gb} GB"
        docker: 'quay.io/biocontainers/ngsderive:3.3.2--pyhdfd78af_0'
        maxRetries: max_retries
    }
}

task read_length {
    input {
        File bam
        File bam_index
        String outfile_name = basename(bam, ".bam") + ".readlength.tsv"
        Float majority_vote_cutoff = 0.7
        Int num_reads = -1
        Int memory_gb = 4
        Int modify_disk_size_gb = 0
        Int max_retries = 1
    }
    
    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        # localize BAM and BAI to CWD
        # some backends prevent writing to the inputs directories
        # to accomodate this, create symlinks in CWD
        CWD_BAM=~{basename(bam)}
        ln -s ~{bam} "$CWD_BAM"
        ln -s ~{bam_index} "$CWD_BAM".bai

        ngsderive readlen --verbose \
            -c ~{majority_vote_cutoff} \
            -n ~{num_reads} \
            "$CWD_BAM" \
            > ~{outfile_name}

        rm "$CWD_BAM" "$CWD_BAM".bai
    >>>

    output {
        File read_length_file = outfile_name
    }

    runtime {
        memory: "~{memory_gb} GB"
        disk: "~{disk_size_gb} GB"
        docker: 'quay.io/biocontainers/ngsderive:3.3.2--pyhdfd78af_0'
        maxRetries: max_retries
    }
}

task encoding {
    input {
        Array[File] ngs_files
        String outfile_name
        Int num_reads = 1000000
        Int memory_gb = 4
        Int modify_disk_size_gb = 0
        Int max_retries = 1
    }

    Float files_size = size(ngs_files, "GiB")
    Int disk_size_gb = ceil(files_size) + 10 + modify_disk_size_gb
 
    command <<<
        set -euo pipefail

        ngsderive encoding --verbose \
            -n ~{num_reads} \
            ~{sep(" ", ngs_files)} \
            > ~{outfile_name}
        
        ENCODING_FILE="~{outfile_name}" python - <<END
import os  # lint-check: ignore

encoding_file = open(os.environ['ENCODING_FILE'], 'r')
encoding_file.readline()  # discard header
permissive_encoding = ""
for line in encoding_file:
    contents = line.strip().split('\t')
    cur_encoding = contents[2]
    if cur_encoding == "Unkown":
        permissive_encoding = cur_encoding
        break
    if cur_encoding == "Sanger/Illumina 1.8" or permissive_encoding == "Sanger/Illumina 1.8":
        permissive_encoding = "Sanger/Illumina 1.8"
    elif cur_encoding == "Solexa/Illumina 1.0" or permissive_encoding == "Solexa/Illumina 1.0":
        permissive_encoding = "Solexa/Illumina 1.0"
    elif cur_encoding == "Illumina 1.3":
        permissive_encoding = "Illumina 1.3"

outfile = open("encoding.txt", "w")
outfile.write(permissive_encoding)
outfile.close()
    
END
    >>>

    output {
        String inferred_encoding = read_string("encoding.txt")
        File encoding_file = outfile_name
    }

    runtime {
        memory: "~{memory_gb} GB"
        disk: "~{disk_size_gb} GB"
        docker: 'quay.io/biocontainers/ngsderive:3.3.2--pyhdfd78af_0'
        maxRetries: max_retries
    }
}

task junction_annotation {
    input {
        File bam
        File bam_index
        File gtf
        String prefix = basename(bam, ".bam")
        Int min_intron = 50
        Int min_mapq = 30
        Int min_reads = 2
        Int fuzzy_junction_match_range = 0
        Int memory_gb = 56  # TODO make this dynamic
        Int modify_disk_size_gb = 0
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        # localize BAM and BAI to CWD
        # some backends prevent writing to the inputs directories
        # to accomodate this, create symlinks in CWD
        CWD_BAM=~{basename(bam)}
        ln -s ~{bam} "$CWD_BAM"
        ln -s ~{bam_index} "$CWD_BAM".bai

        ngsderive junction-annotation --verbose \
            -g ~{gtf} \
            -i ~{min_intron} \
            -q ~{min_mapq} \
            -m ~{min_reads} \
            -k ~{fuzzy_junction_match_range} \
            -o ~{prefix}.junction_summary.tsv \
            "$CWD_BAM"

        # junction-annotation accepts multiple BAMs, and allows for
        # renaming the cohort level summary report, but not the BAM
        # level junctions file. So we rename it here to match with prefix.
        mv "$(basename ~{bam}.junctions.tsv)" "~{prefix}.junctions.tsv"
        gzip ~{prefix}.junctions.tsv

        rm "$CWD_BAM" "$CWD_BAM".bai
    >>>

    output {
        File junction_summary = "~{prefix}.junction_summary.tsv"
        File junctions = "~{prefix}.junctions.tsv.gz"
    }

    runtime {
        memory: "~{memory_gb} GB"
        disk: "~{disk_size_gb} GB"
        docker: 'quay.io/biocontainers/ngsderive:3.3.2--pyhdfd78af_0'
        maxRetries: max_retries
    }
}

task endedness {
    meta {
        description: "This WDL task derives the endedness of the input BAM file"
    }

    parameter_meta {
        bam: "Input BAM format file to derive endedness from"
        outfile_name: "Name for the endedness TSV file"
        lenient: "Return a zero exit code on unknown results"
        calc_rpt: "Calculate and output Reads-Per-Template. This will produce a more sophisticated estimate for endedness, but uses substantially more memory (can reach up to 60-70% of BAM size in memory consumption)."
        round_rpt: "Round RPT to the nearest INT before comparing to expected values. Appropriate if using `--num-reads` > 0."
        split_by_rg: "Contain one entry per read group"
        paired_deviance: "Distance from 0.5 split between number of f+l- reads and f-l+ reads allowed to be called 'Paired-End'. Default of `0.0` only appropriate if the whole file is being processed."
        num_reads: "How many reads to analyze from the start of the file. Any n < 1 to parse whole file."
        modify_memory_gb: "Add to or subtract from dynamic memory allocation. Default memory is determined by value of `calc_rpt` and the size of the input. Specified in GB."
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
        max_retries: "Number of times to retry in case of failure"
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
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int memory_gb = if calc_rpt then (
        ceil(bam_size * 1.75) + 4 + modify_memory_gb
    ) else 4
    Int disk_size_gb = ceil(bam_size) + 10 + modify_disk_size_gb

    command <<<
        ngsderive endedness --verbose \
            ~{if lenient then "--lenient" else ""} \
            ~{if calc_rpt then "-r" else ""} \
            ~{if round_rpt then "--round-rpt" else ""} \
            ~{if split_by_rg then "--split-by-rg" else ""} \
            --paired-deviance ~{paired_deviance} \
            -n ~{num_reads} \
            ~{bam} \
            > ~{outfile_name}
    >>>

    output {
        File endedness_file = outfile_name
    }

    runtime {
        memory: "~{memory_gb} GB"
        disk: "~{disk_size_gb} GB"
        docker: 'quay.io/biocontainers/ngsderive:3.3.2--pyhdfd78af_0'
        maxRetries: max_retries
    }
}
