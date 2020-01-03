## Description: 
##
## This WDL tool wraps the PicardTools library (https://broadinstitute.github.io/picard/).
## PicardTools is a set of Java tools for manipulating sequencing data. 

task mark_duplicates {
    File bam
    String prefix = basename(bam, ".bam")

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
        memory: "50G"
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }

    output {
        File out = "${prefix}.duplicates.bam"
    }
    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
        author: "Andrew Frantz"
        email: "andrew.frantz@stjude.org"
        description: "This WDL tool marks duplicate reads in the input BAM file using Picard."
    }
    parameter_meta {
        bam: "Input BAM format file to generate coverage for"
    }
}

task validate_bam {
    File bam
    
    command {
        picard ValidateSamFile I=${bam} \
            IGNORE=INVALID_PLATFORM_VALUE > stdout.txt
    }

    runtime {
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }

    output {
        String out = read_string("stdout.txt")
    }
    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
        author: "Andrew Frantz"
        email: "andrew.frantz@stjude.org"
        description: "This WDL tool validates the input BAM file for correct formatting using Picard."
    }
    parameter_meta {
        bam: "Input BAM format file to generate coverage for"
    }
}

task bam_to_fastq {
    File bam
    String prefix = basename(bam, ".bam")

    command {
        picard SamToFastq INPUT=${bam} \
            FASTQ=${prefix}_R1.fastq \
            SECOND_END_FASTQ=${prefix}_R2.fastq \
            RE_REVERSE=true
    }

    runtime{
        memory: "25G"	
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }

    output {
        File read1 = "${prefix}_R1.fastq"
        File read2 = "${prefix}_R2.fastq"
    }
    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
        author: "Andrew Frantz"
        email: "andrew.frantz@stjude.org"
        description: "This WDL tool converts the input BAM file into paired FastQ format files."
    }
    parameter_meta {
        bam: "Input BAM format file to generate coverage for"
    }
}
