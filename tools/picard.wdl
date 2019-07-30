task mark_duplicates {
    File bam
    String basename = basename(bam, ".bam")

    command {
        picard MarkDuplicates I=${bam} \
            O=${basename}.duplicates.bam \
            VALIDATION_STRINGENCY=SILENT \
            CREATE_INDEX=false \
            CREATE_MD5_FILE=false \
            COMPRESSION_LEVEL=5 \
            METRICS_FILE=${basename}.metrics.txt
    }

    output {
        File out = "${basename}.duplicates.bam"
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
    String basename = basename(bam, ".bam")

    command {
        picard SamToFastq INPUT=${bam} \
            FASTQ=${basename}_R1.fastq \
            SECOND_END_FASTQ=${basename}_R2.fastq \
            RE_REVERSE=true
    }

    output {
        File read1 = "${basename}_R1.fastq"
        File read2 = "${basename}_R2.fastq"
    }
}
