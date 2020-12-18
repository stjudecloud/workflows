## # BWA
##
## This WDL tool wraps [BWA](https://github.com/lh3/bwa).
## BWA aligns sequencing reads to a reference genome.

version 1.0

task bwa_aln {
    input {
        File fastq
        String output_bam = basename(fastq, ".fq.gz") + ".bam"
        File bwadb_tar_gz
        Int ncpu = 1
        Int memory_gb = 5
        Int? disk_size_gb
        Int max_retries = 1
        Boolean detect_nproc = false
    }

    String parsed_detect_nproc = if detect_nproc then "true" else ""
    Float input_fastq_size = size(fastq, "GiB")
    Float reference_size = size(bwadb_tar_gz, "GiB")
    Int disk_size = select_first([disk_size_gb, ceil((input_fastq_size * 2) + (reference_size * 2))])

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if [ -n ~{parsed_detect_nproc} ]
        then
            n_cores=$(nproc)
        fi

        mkdir bwa
        tar -C bwa -xzf ~{bwadb_tar_gz}
        PREFIX=$(basename bwa/*.ann ".ann")

        bwa aln bwa/$PREFIX ~{fastq} > sai

        bwa samse bwa/$PREFIX sai ~{fastq} | java.sh org.stjude.compbio.sam.TweakSam -V SILENT -G 4 -o ~{output_bam}        
    >>>

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        cpu: ncpu
        docker: 'stjudecloud/bwa:1.0.0'
        maxRetries: max_retries
    }

    output {
        File bam = output_bam
    }

    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
        description: "This WDL tool converts a SAM file to FastQ format." 
    }

    parameter_meta {
        fastq: "Input FastQ file to align with bwa"
    }
}


task bwa_mem {
    input {
        File fastq
        String output_bam = basename(fastq, ".fq.gz") + ".bam"
        File bwadb_tar_gz
        Int ncpu = 1
        Int memory_gb = 5
        Int? disk_size_gb
        Int max_retries = 1
        Boolean detect_nproc = false
    }

    String parsed_detect_nproc = if detect_nproc then "true" else ""
    Float input_fastq_size = size(fastq, "GiB")
    Float reference_size = size(bwadb_tar_gz, "GiB")
    Int disk_size = select_first([disk_size_gb, ceil((input_fastq_size * 2) + reference_size)])

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if [ -n ~{parsed_detect_nproc} ]
        then
            n_cores=$(nproc)
        fi

        mkdir bwa
        tar -C bwa -xzf ~{bwadb_tar_gz}
        PREFIX=$(basename bwa/*.ann ".ann")

        bwa mem -t $n_cores bwa/$PREFIX ~{fastq} > ~{output_bam}
    >>>

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        cpu: ncpu
        docker: 'stjudecloud/bwa:1.0.0'
        maxRetries: max_retries
    }

    output {
        File bam = output_bam
    }

    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
        description: "This WDL tool converts a SAM file to FastQ format." 
    }

    parameter_meta {
        fastq: "Input FastQ file to align with bwa"
    }
}