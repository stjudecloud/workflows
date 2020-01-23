## Description:
##
## This WDL tool wraps the PicardTools library (https://broadinstitute.github.io/picard/).
## PicardTools is a set of Java tools for manipulating sequencing data.

version 1.0

task mark_duplicates {
    input {
        File bam
        String prefix = basename(bam, ".bam")
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil((bam_size * 2) + 10)

    command {
        picard MarkDuplicates I=${bam} \
            O=${prefix}.duplicates.bam \
            VALIDATION_STRINGENCY=SILENT \
            CREATE_INDEX=false \
            CREATE_MD5_FILE=false \
            COMPRESSION_LEVEL=5 \
            METRICS_FILE=${prefix}.metrics.txt
    }

    runtime {
        memory: "50 GB"
        disk: disk_size + " GB"
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }

    output {
        File out = "${prefix}.duplicates.bam"
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "This WDL tool marks duplicate reads in the input BAM file using Picard."
    }

    parameter_meta {
        bam: "Input BAM format file to mark duplicates"
    }
}

task validate_bam {
    input {
        File bam
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil((bam_size * 2) + 10)
    
    command {
        picard ValidateSamFile I=${bam} \
            IGNORE=INVALID_PLATFORM_VALUE > stdout.txt
    }

    runtime {
        disk: disk_size + " GB"
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }

    output {
        String out = read_string("stdout.txt")
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "This WDL tool validates the input BAM file for correct formatting using Picard."
    }

    parameter_meta {
        bam: "Input BAM format file to validate"
    }
}

task bam_to_fastq {
    input {
        File bam
        String prefix = basename(bam, ".bam")
        Int memory_gb = 25
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil((bam_size * 4) + 10)

    command {
        picard SamToFastq INPUT=${bam} \
            FASTQ=${prefix}_R1.fastq \
            SECOND_END_FASTQ=${prefix}_R2.fastq \
            RE_REVERSE=true \
            VALIDATION_STRINGENCY=SILENT
    }

    runtime{
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }

    output {
        File read1 = "${prefix}_R1.fastq"
        File read2 = "${prefix}_R2.fastq"
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "This WDL tool converts the input BAM file into paired FastQ format files."
    }

    parameter_meta {
        bam: "Input BAM format file to convert to FastQ"
    }
}

task sort {
    input {
        File bam
        String sort_order = "coordinate"
        String output_filename = basename(bam, ".bam") + ".sorted.bam"
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil((bam_size * 4) + 10)

    command {
        picard SortSam I=${bam} \
           O=${output_filename} \
           SO=${sort_order} 
    }
    runtime {
        memory: "25 GB"
        disk: disk_size + " GB"
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }
    output {
        File sorted_bam = output_filename
    }
    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
        description: "This WDL tool sorts the input BAM file."
    }
    parameter_meta {
        bam: "Input BAM format file to sort"
    }
}
