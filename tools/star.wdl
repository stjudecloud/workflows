## # STAR
##
## This WDL file wraps the [STAR aligner](https://github.com/alexdobin/STAR).
## STAR is an RNA-Seq aligner.

version 1.1

task build_star_db {
    meta {
        description: "This WDL task runs STAR's build command to generate a STAR format reference for alignment." 
    }

    parameter_meta {
        reference_fasta: "The FASTA format reference file for the genome"
        gtf: "GTF format feature file"
    }

    input {
        File reference_fasta
        File gtf
        String db_name = "star_db"
        Boolean use_all_cores = false
        Int ncpu = 1
        Int memory_gb = 50
        Int modify_disk_size_gb = 0
        Int max_retries = 1
    }

    String star_db_tar_gz = db_name + ".tar.gz"

    Float reference_fasta_size = size(reference_fasta, "GiB")
    Float gtf_size = size(gtf, "GiB")
    Int disk_size_gb = (
        ceil((reference_fasta_size + gtf_size) * 3) + 10 + modify_disk_size_gb
    )

    # Leave 2GB as system overhead
    String memory_limit_bytes = "~{memory_gb - 2}000000000"

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi

        gtf_name=~{basename(gtf, ".gz")}
        gunzip -c ~{gtf} > "$gtf_name" || ln -s ~{gtf} "$gtf_name"

        ref_fasta=~{basename(reference_fasta, ".gz")}
        gunzip -c ~{reference_fasta} > "$ref_fasta" \
            || ln -s ~{reference_fasta} "$ref_fasta"
        
        mkdir ~{db_name};
        STAR --runMode genomeGenerate \
            --genomeDir ~{db_name} \
            --runThreadN "$n_cores" \
            --limitGenomeGenerateRAM ~{memory_limit_bytes} \
            --genomeFastaFiles "$ref_fasta" \
            --sjdbGTFfile "$gtf_name" \
            --sjdbOverhang 125

        rm "$gtf_name" "$ref_fasta"

        tar -czf ~{star_db_tar_gz} ~{db_name}
    >>>

    output {
        File star_db = star_db_tar_gz
    }

    runtime {
        cpu: ncpu
        memory: "~{memory_gb} GB"
        disk: "~{disk_size_gb} GB"
        docker: 'ghcr.io/stjudecloud/star:2.7.10a-0'
        maxRetries: max_retries
    }
}

task alignment {
    meta {
        description: "This WDL task runs the STAR aligner on a set of RNA-Seq FastQ files."
    }

    parameter_meta {
        read_one_fastqs: "An array of FastQ files containing read one information"
        star_db_tar_gz: "A gzipped TAR file containing the STAR reference files. The name of the root directory which was archived must match the archive's filename without the `.tar.gz` extension."
        read_groups: "A string containing the read group information to output in the BAM file. If including multiple read group fields per-read group, they should be space delimited. Read groups should be comma separated, with a space on each side (e.g. ' , '). The ID field must come first for each read group and must match the basename of a fastq file (up to the first period). Example: `ID:rg1 PU:flowcell1.lane1 SM:sample1 PL:illumina LB:sample1_lib1 , ID:rg2 PU:flowcell1.lane2 SM:sample1 PL:illumina LB:sample1_lib1`"
        read_two_fastqs: "An array of FastQ files containing read two information"
    }

    input {
        Array[File] read_one_fastqs
        File star_db_tar_gz
        String prefix
        String? read_groups
        Array[File] read_two_fastqs = []
        Boolean use_all_cores = false
        Int ncpu = 1
        Int memory_gb = 50
        Int modify_disk_size_gb = 0
        Int max_retries = 1
    }
    
    String star_db_dir = basename(star_db_tar_gz, ".tar.gz")

    # Leave 2GB as system overhead
    String memory_limit_bytes = "~{memory_gb - 2}000000000"

    Float read_one_fastqs_size = size(read_one_fastqs, "GiB")
    Float read_two_fastqs_size = size(read_two_fastqs, "GiB")
    Float star_db_tar_gz_size = size(star_db_tar_gz, "GiB")
    Int disk_size_gb = (
        (
            ceil(read_one_fastqs_size + read_two_fastqs_size + star_db_tar_gz_size) * 3
        ) + 10 + modify_disk_size_gb
    )

    Array[File] empty_array = []  # odd construction forced by WDL v1.0 spec

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi

        tar -xzf ~{star_db_tar_gz}

        # TODO rework `sort_star_input.py` to avoid this spaghetti logic
        if [ -n "~{if defined(read_groups) then read_groups else ""}" ]
        then
            if [ -n "~{if read_two_fastqs != empty_array then "read_two_fastqs" else ""}" ]
            then
                python3 /home/sort_star_input.py \
                    --read_one_fastqs "~{sep(',', read_one_fastqs)}" \
                    --read_two_fastqs "~{sep(',', read_two_fastqs)}" \
                    --read_groups "~{read_groups}"
            else
                python3 /home/sort_star_input.py \
                    --read_one_fastqs "~{sep(',', read_one_fastqs)}" \
                    --read_groups "~{read_groups}"
            fi
        else 
            if [ -n "~{if read_two_fastqs != empty_array then "read_two_fastqs" else ""}" ]
            then
                python3 /home/sort_star_input.py \
                    --read_one_fastqs "~{sep(',', read_one_fastqs)}" \
                    --read_two_fastqs "~{sep(',', read_two_fastqs)}"
            else
                python3 /home/sort_star_input.py \
                    --read_one_fastqs "~{sep(',', read_one_fastqs)}"
            fi
        fi

        STAR --readFilesIn \
                $(cat read_one_fastqs_sorted.txt) \
                $(cat read_two_fastqs_sorted.txt) \
             --readFilesCommand "gunzip -c" \
             --genomeDir ~{star_db_dir} \
             --runThreadN "$n_cores" \
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
             --outFileNamePrefix ~{prefix + "."} \
             --twopassMode Basic \
             --limitBAMsortRAM ~{memory_limit_bytes} \
             --outSAMattrRGline $(cat read_groups_sorted.txt)
    >>>

    output {
        File star_log = prefix + ".Log.final.out"
        File star_bam = prefix + ".Aligned.out.bam"
    }

    runtime {
        cpu: ncpu
        memory: "~{memory_gb} GB"
        disk: "~{disk_size_gb} GB"
        docker: 'ghcr.io/stjudecloud/star:2.7.10a-0'
        maxRetries: max_retries
    }
}
