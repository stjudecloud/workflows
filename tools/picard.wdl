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
            # METRICS_FILE=$METRICS_FILE
    }

    output {
        File out = "${basename}.duplicates.bam"
    }
}

task validate_bam {
    File bam
    
    command {
        picard ValidateSamFile I=${bam} \
            IGNORE=INVALID_PLATFORM_VALUE
    }
}
