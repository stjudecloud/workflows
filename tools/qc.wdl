task quickcheck {
    File bam

    command {
        samtools quickcheck ${bam}
    }
}

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

task fastqc {
    File bam
    Int ncpu
    String basename = basename(bam, ".bam")

    command {
        fastqc -f bam \
            -o ${basename}_fastqc_results \
            -t ${ncpu} \
            ${bam}; \
        cd ${basename}_fastqc_results; \
        dir=`pwd`; \
        ln *.gz "$dir"; \
        ls > file-list.txt
    }

    output {
        Array[File] files = read_lines("file-list.txt")
    }
}

task qualimap {
    File bam
    Int ncpu
    String basename = basename(bam, ".bam")

    command {
        qualimap bamqc -bam ${bam} \
            -outdir ${basename}_qualimap_results \
            -nt ${ncpu}; \
        cd ${basename}_fastqc_results; \
        dir=`pwd`; \
        ln *.gz "$dir"; \
        ls > file-list.txt
    }

    output {
        Array[File] files = read_lines("file-list.txt")
    }
}
