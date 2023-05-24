## # ngsderive
##
## This WDL file wraps the [ngsderive package](https://github.com/stjudecloud/ngsderive)

version 1.0

task infer_strandedness {
    input {
        File bam
        File bam_index
        File gtf
        Int min_reads_per_gene = 10
        Int num_genes = 1000
        Int min_mapq = 30
        Int memory_gb = 5
        Int max_retries = 1
    }

    String out_file = basename(bam, ".bam") + ".strandedness.txt"
    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil((bam_size) + 10)
 
    command {
        set -euo pipefail

        if [ ! -e ~{bam}.bai ]
        then 
            ln -s ~{bam_index} ~{bam}.bai
        fi

        ngsderive strandedness --verbose \
            -m ~{min_reads_per_gene} \
            -n ~{num_genes} \
            -q ~{min_mapq} \
            ~{bam} \
            -g ~{gtf} \
            > ~{out_file}
        awk 'NR > 1' ~{out_file} | cut -d$'\t' -f5 > strandedness.txt
    }

    runtime {
        disk: disk_size + " GB"
        docker: 'quay.io/biocontainers/ngsderive:2.4.0--pyhdfd78af_0'
        memory: memory_gb + " GB"
        maxRetries: max_retries
    }

    output {
        String strandedness = read_string("strandedness.txt")
        File strandedness_file = out_file
    }
}

task instrument {
    input {
        File bam
        Int num_samples = 10000
        Int max_retries = 1
    }

    String out_file = basename(bam, ".bam") + ".instrument.txt"
    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil((bam_size) + 10)

    command {
        ngsderive instrument --verbose \
            -n ~{num_samples} \
            ~{bam} \
            > ~{out_file}
    }

    runtime {
        memory: "4 GB"
        disk: disk_size + " GB"
        docker: 'quay.io/biocontainers/ngsderive:2.4.0--pyhdfd78af_0'
        maxRetries: max_retries
    }

    output {
        File instrument_file = out_file
    }
}

task read_length {
    input {
        File bam
        File bam_index
        Float majority_vote_cutoff = 0.7
        Int num_samples = 10000
        Int memory_gb = 5
        Int max_retries = 1
    }

    String out_file = basename(bam, ".bam") + ".readlength.txt"
    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil(bam_size + 10)

    command {
        set -euo pipefail

        if [ ! -e ~{bam}.bai ]
        then 
            ln -s ~{bam_index} ~{bam}.bai
        fi

        ngsderive readlen --verbose \
            -c ~{majority_vote_cutoff} \
            -n ~{num_samples} \
            ~{bam} \
            > ~{out_file}
    }

    runtime {
        disk: disk_size + " GB"
        memory: memory_gb + " GB"
        docker: 'quay.io/biocontainers/ngsderive:2.4.0--pyhdfd78af_0'
        maxRetries: max_retries
    }

    output {
        File read_length_file = out_file
    }

}

task encoding {
    input {
        Array[File] ngs_files
        String prefix
        Int num_samples = 1000000
        Int memory_gb = 5
        Int max_retries = 1
    }

    String out_file = prefix + ".encoding.txt"
    Float files_size = size(ngs_files, "GiB")
    Int disk_size = ceil((files_size) + 5)
 
    command <<<
        set -euo pipefail

        ngsderive encoding --verbose \
            -n ~{num_samples} \
            ~{sep=' ' ngs_files} \
            > ~{out_file}
        
        ENCODING_FILE="~{out_file}" python - <<END
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

    runtime {
        disk: disk_size + " GB"
        memory: memory_gb + " GB"
        docker: 'quay.io/biocontainers/ngsderive:2.4.0--pyhdfd78af_0'
        maxRetries: max_retries
    }

    output {
        String inferred_encoding = read_string("encoding.txt")
        File encoding_file = out_file
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
        Int memory_gb = 35
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil((bam_size) + 10)

    command {
        set -euo pipefail

        if [ ! -e ~{bam}.bai ]
        then 
            ln -s ~{bam_index} ~{bam}.bai
        fi

        ngsderive junction-annotation --verbose \
            -g ~{gtf} \
            -i ~{min_intron} \
            -q ~{min_mapq} \
            -m ~{min_reads} \
            -k ~{fuzzy_junction_match_range} \
            -o ~{prefix}.junction_summary.txt \
            ~{bam}

        mv "$(basename ~{bam}.junctions.tsv)" "~{prefix}.junctions.tsv"
        gzip ~{prefix}.junctions.tsv
    }

    runtime {
        disk: disk_size + " GB"
        memory: memory_gb + " GB"
        docker: 'quay.io/biocontainers/ngsderive:2.4.0--pyhdfd78af_0'
        maxRetries: max_retries
    }

    output {
        File junction_summary = "~{prefix}.junction_summary.txt"
        File junctions = "~{prefix}.junctions.tsv.gz"
    }
}
