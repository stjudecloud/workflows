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
    }

    output {
        File out = "${prefix}.duplicates.bam"
    }
}

task validate_bam {
    File bam
    
    command {
        picard ValidateSamFile I=${bam} \
            IGNORE=INVALID_PLATFORM_VALUE > stdout.txt
    }

    output {
        String out = read_string("stdout.txt")
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
    }

    output {
        File read1 = "${prefix}_R1.fastq"
        File read2 = "${prefix}_R2.fastq"
    }
}
