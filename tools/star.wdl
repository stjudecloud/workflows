task print_version {
    command {
        STAR --version
    }

    output {
        String out = read_string(stdout())
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