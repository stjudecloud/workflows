## # STAR
##
## This WDL tool wraps the [STAR aligner](https://github.com/alexdobin/STAR).
## STAR is an RNA-seq aligner.

version 1.0

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
        Boolean detect_nproc = false
    }

    String parsed_detect_nproc = if detect_nproc then "true" else ""
    String stardb_out_name = stardb_dir_name + ".tar.gz"
    Float reference_fasta_size = size(reference_fasta, "GiB")
    Float gencode_gtf_size = size(gencode_gtf, "GiB")
    Int disk_size = select_first([disk_size_gb, ceil(((reference_fasta_size + gencode_gtf_size) * 3) + 10)])

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if [ -n ~{parsed_detect_nproc} ]
        then
            n_cores=$(nproc)
        fi

        orig_gtf=~{gencode_gtf}
        gtf=$(basename "${orig_gtf%.gz}")
        gunzip -c ~{gencode_gtf} > $gtf || cp ~{gencode_gtf} $gtf

        orig_fasta=~{reference_fasta}
        ref_fasta=$(basename "${orig_fasta%.gz}")
        gunzip -c ~{reference_fasta} > $ref_fasta || cp ~{reference_fasta} $ref_fasta
        
        mkdir ~{stardb_dir_name};
        STAR --runMode genomeGenerate \
            --genomeDir ~{stardb_dir_name} \
            --runThreadN $n_cores \
            --limitGenomeGenerateRAM ~{ram_limit} \
            --genomeFastaFiles $ref_fasta \
            --sjdbGTFfile $gtf \
            --sjdbOverhang 125
        rm $gtf $ref_fasta
        tar -czf ~{stardb_out_name} ~{stardb_dir_name}
    >>>

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        cpu: ncpu
        docker: 'stjudecloud/star:1.1.1'
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
        Array[File]? read_two_fastqs
        File stardb_tar_gz
        String output_prefix
        String? read_groups
        Int memory_gb = 50
        Int? disk_size_gb
        Int max_retries = 1
        Boolean detect_nproc = false
    }
    
    String parsed_detect_nproc = if detect_nproc then "true" else ""
    String stardb_dir = basename(stardb_tar_gz, ".tar.gz")
    Float read_one_fastqs_size = size(read_one_fastqs, "GiB")
    Float read_two_fastqs_size = size(select_first([read_two_fastqs, []]), "GiB")
    Float stardb_tar_gz_size = size(stardb_tar_gz, "GiB")
    Int disk_size = select_first([disk_size_gb, ceil(((read_one_fastqs_size + read_two_fastqs_size + stardb_tar_gz_size) * 3) + 10)])

    command {
        set -euo pipefail

        n_cores=~{ncpu}
        if [ -n ~{parsed_detect_nproc} ]
        then
            n_cores=$(nproc)
        fi

        tar -xzf ~{stardb_tar_gz};

        if [ -n "~{if defined(read_groups) then read_groups else ""}" ]
        then
            if [ -n "~{if defined(read_two_fastqs) then "read_two_fastqs" else ""}" ]
            then
                python /home/sort_star_input.py \
                    --read_one_fastqs "~{sep=',' read_one_fastqs}" \
                    --read_two_fastqs "~{sep=',' read_two_fastqs}" \
                    --read_groups "~{read_groups}"
            else
                python /home/sort_star_input.py \
                    --read_one_fastqs "~{sep=',' read_one_fastqs}" \
                    --read_groups "~{read_groups}"
            fi
        else 
            if [ -n "~{if defined(read_two_fastqs) then "read_two_fastqs" else ""}" ]
            then
                python /home/sort_star_input.py \
                    --read_one_fastqs "~{sep=',' read_one_fastqs}" \
                    --read_two_fastqs "~{sep=',' read_two_fastqs}"
            else
                python /home/sort_star_input.py \
                    --read_one_fastqs "~{sep=',' read_one_fastqs}"
            fi
        fi

        STAR --readFilesIn $(cat read_one_fastqs_sorted.txt) $(cat read_two_fastqs_sorted.txt) \
             --readFilesCommand "gunzip -c" \
             --genomeDir ~{stardb_dir} \
             --runThreadN $n_cores \
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
        docker: 'stjudecloud/star:1.1.1'
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
        read_groups: "A string containing the read group information to output in the BAM file. If including multiple read group fields per-read group, they should be space delimited. Read groups should be comma separated, with a space on each side (e.g. ' , '). The ID field must come first for each read group and must match the basename of a fastq file (up to the first period). Example: `ID:rg1 PU:flowcell1.lane1 SM:sample1 PL:illumina LB:sample1_lib1 , ID:rg2 PU:flowcell1.lane2 SM:sample1 PL:illumina LB:sample1_lib1`"
    }
}
