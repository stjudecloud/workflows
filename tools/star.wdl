## Description:
##
## This WDL tool wraps the STAR aligner (https://github.com/alexdobin/STAR).
## STAR is an RNA-seq aligner.

version 1.0

task print_version {
    command {
        STAR --version
    }

    runtime {
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }

    output {
        String out = read_string(stdout())
    }
}

task build_db {
    input {
        Int ncpu
        File reference_fasta
        File gencode_gtf
        String stardb_dir_name
        String stardb_zip_name = stardb_dir_name + ".zip"
        String ram_limit = "45000000000"
    }
    Float reference_fasta_size = size(reference_fasta, "GiB")
    Float gencode_gtf_size = size(gencode_gtf, "GiB")
    Int disk_size = ceil(((reference_fasta_size + gencode_gtf_size) * 3) + 10)

    command {
        mkdir ${stardb_dir_name};
        STAR --runMode genomeGenerate \
            --genomeDir ${stardb_dir_name} \
            --runThreadN ${ncpu} \
            --limitGenomeGenerateRAM ${ram_limit} \
            --genomeFastaFiles ${reference_fasta} \
            --sjdbGTFfile ${gencode_gtf} \
            --sjdbOverhang 125
        zip -r ${stardb_zip_name} ${stardb_dir_name}
    }

    runtime {
        memory: "50 GB"
        disk: disk_size + " GB"
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }

    output {
        File zip = stardb_zip_name
    }
}

task alignment {
    input {
      Array[File] read_one_fastqs
      Array[File] read_two_fastqs
      File stardb_zip
      String stardb_dir = basename(stardb_zip, ".zip")
      String output_prefix
      String? read_groups
    }
    Float read_one_fastqs_size = size(read_one_fastqs, "GiB")
    Float read_two_fastqs_size = size(read_two_fastqs, "GiB")
    Float stardb_zip_size = size(stardb_zip, "GiB")
    Int disk_size = ceil(((read_one_fastqs_size + read_two_fastqs_size + stardb_zip_size) * 3) + 10)

    command {
        unzip ${stardb_zip};
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
             --twopassMode Basic \
             --limitBAMsortRAM 58000000000 \
             ${"--outSAMattrRGline " + read_groups}
    }

    runtime {
        memory: "50 GB"
        disk: disk_size + " GB"
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }

    output {
       File star_bam = output_prefix + "Aligned.sortedByCoord.out.bam"
    }
}
