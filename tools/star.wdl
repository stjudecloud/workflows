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
    String stardb_dir_name

    command {
        mkdir ${stardb_dir_name};
        STAR --runMode genomeGenerate \
            --genomeDir ${stardb_dir_name} \
            --runThreadN ${ncpu} \
            --genomeFastaFiles ${reference_fasta} \
            --sjdbGTFfile ${gencode_gtf} \
            --sjdbOverhang 125
    }

    output {
        File out_dir = stardb_dir_name
    }
}

task alignment {
    Array[File] read_one_fastqs
    Array[File] read_two_fastqs
    File stardb_dir
    String output_prefix

    command {
        STAR --readFilesIn ${sep=',' read_one_fastqs} ${sep=',' read_two_fastqs} \
             --genomeDir ${stardb_dir} \
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
             --outFileNamePrefix ${output_prefix} \
             --twopassMode Basic
    }
}