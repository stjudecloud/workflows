task print_version {
    command {
        STAR --version
    }

    output {
        String out = read_string(stdout())
    }
}

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

task build_db {
    Int ncpu
    File reference_fasta
    File gencode_gtf
    File stardb

    command {
        STAR --runMode genomeGenerate \
            --genomeDir ${stardb} \
            --runThreadN ${ncpu} \
            --genomeFastaFiles ${reference_fasta} \
            --sjdbGTFfile ${gencode_gtf} \
            --sjdbOverhang 125
    }
}

task alignment {
    File read_one_fastq
    File read_two_fastq
    File stardb
    String outprefix

    command {
        echo STAR --readFilesIn ${read_one_fastq} ${read_two_fastq} \
             --genomeDir ${stardb} \
             --outSAMunmapped Within \
             --outSAMstrandField intronMotif \
             --outSAMtype BAM SortedByCoordinate \
             --outSAMattributes NH HI AS nM NM MD XS \
             --outFilterMultimapScoreRange 1 \
             --outFilterMultimapNmax 20 \
             --outFilterMismatchNmax 10 \
             --alignIntronMax 500000 \
             --alignMatesGapMax 1000000 \
             --sjdbScore 2 \
             --alignSJDBoverhangMin 1 \
             --outFilterMatchNminOverLread 0.66 \
             --outFilterScoreMinOverLread 0.66 \
             --outFileNamePrefix ${outprefix} \
             --twopassMode Basic
    }
}