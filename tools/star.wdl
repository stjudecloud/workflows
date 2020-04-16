## # STAR
##
## This WDL tool wraps the [STAR aligner](https://github.com/alexdobin/STAR).
## STAR is an RNA-seq aligner.

version 1.0

task star_print_version {
    command {
        STAR --version
    }

    runtime {
        docker: 'stjudecloud/star:1.0.0'
    }

    output {
        String out = read_string(stdout())
    }
}

task build_db {
    input {
        Int ncpu = 1
        File reference_fasta
        File gencode_gtf
        String stardb_dir_name
        String ram_limit = "45000000000" # This value is too large to be an Int type, so we store it as a string
        Int memory_gb = 50
        Int? disk_size_gb
        Int max_retries = 1
    }

    String stardb_out_name = stardb_dir_name + ".tar.gz"
    Float reference_fasta_size = size(reference_fasta, "GiB")
    Float gencode_gtf_size = size(gencode_gtf, "GiB")
    Int disk_size = select_first([disk_size_gb, ceil(((reference_fasta_size + gencode_gtf_size) * 3) + 10)])

    command {
        mkdir ~{stardb_dir_name};
        STAR --runMode genomeGenerate \
            --genomeDir ~{stardb_dir_name} \
            --runThreadN ~{ncpu} \
            --limitGenomeGenerateRAM ~{ram_limit} \
            --genomeFastaFiles ~{reference_fasta} \
            --sjdbGTFfile ~{gencode_gtf} \
            --sjdbOverhang 125
        tar -czf ~{stardb_out_name} ~{stardb_dir_name}
    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        cpu: ncpu
        docker: 'stjudecloud/star:1.0.0'
        maxRetries: max_retries
    }

    output {
        File stardb_out = stardb_out_name
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "This WDL tool runs STAR's build command to generate a STAR format reference for alignment." 
    }

    parameter_meta {
        reference_fasta: "The FASTA format reference file for the genome"
        gencode_gtf: "GTF format feature file with Gencode features"
    }
}

task alignment {
    input {
        Int ncpu = 1
        Array[File] read_one_fastqs
        Array[File] read_two_fastqs
        File stardb_tar_gz
        String output_prefix
        String? read_groups
        Int memory_gb = 50
        Int? disk_size_gb
        Int max_retries = 1
    }
    
    String stardb_dir = basename(stardb_tar_gz, ".tar.gz")
    Float read_one_fastqs_size = size(read_one_fastqs, "GiB")
    Float read_two_fastqs_size = size(read_two_fastqs, "GiB")
    Float stardb_tar_gz_size = size(stardb_tar_gz, "GiB")
    Int disk_size = select_first([disk_size_gb, ceil(((read_one_fastqs_size + read_two_fastqs_size + stardb_tar_gz_size) * 3) + 10)])

    command {
        tar -xzf ~{stardb_tar_gz};

        python /home/sort_star_input.py \
            --read_one_fastqs "~{sep=',' read_one_fastqs}" \
            --read_two_fastqs "~{sep=',' read_two_fastqs}" \
            --read_groups "~{read_groups}"

        STAR --readFilesIn $(cat read_one_fastqs_sorted.txt) $(cat read_two_fastqs_sorted.txt) \
             --genomeDir ~{stardb_dir} \
             --runThreadN ~{ncpu} \
             --outSAMunmapped Within \
             --outSAMstrandField intronMotif \
             --outSAMtype BAM Unsorted \
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
             --outFileNamePrefix ~{output_prefix + "."} \
             --twopassMode Basic \
             --limitBAMsortRAM ~{(memory_gb - 2) + "000000000"} \
             --outSAMattrRGline $(cat read_groups_sorted.txt)
    }

    runtime {
        cpu: ncpu
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        docker: 'stjudecloud/star:1.0.0'
        maxRetries: max_retries
    }

    output {
        File star_log = output_prefix + ".Log.final.out"
        File star_bam = output_prefix + ".Aligned.out.bam"
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "This WDL tool generates a FastQC quality control metrics report for the input BAM file."
    }

    parameter_meta {
        read_one_fastqs: "An array of FastQ files containing read one information"
        read_two_fastqs: "An array of FastQ files containing read two information in the same order as the read one FastQ"
        stardb_tar_gz: "A gzipped TAR file containing the STAR reference files"
        read_groups: "A string containing the read group information to output in the BAM file"
    }
}
