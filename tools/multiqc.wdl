## # MultiQC
##
## This WDL tool wraps the [MultiQC](https://multiqc.info/) tool.
## MultiQC aggregates quality control results for bioinformatics.

version 1.0

task multiqc {
    input {
        File sorted_bam
        File validate_sam_file
        File qualimap_bamqc
        File? qualimap_rnaseq
        File fastqc
        File? fastq_screen
        File flagstat_file
        File? star_log
        Int max_retries = 1
        Int memory_gb = 5
    }

    String out_directory = basename(sorted_bam, ".bam")
    String out_tar_gz = out_directory + ".tar.gz"
    Float star_size = size(sorted_bam, "GiB")
    Int disk_size = ceil((star_size * 4) + 10)

    String rnaseq = if defined(qualimap_rnaseq) then "true" else ""
    String fastq_screen_tar = select_first([fastq_screen, []])

    command {
        set -eo pipefail
        
        # set ENV variables for `multiqc`
        export LC_ALL=C.UTF-8
        export LANG=C.UTF-8
        
        echo ~{sorted_bam} > file_list.txt
        echo ~{validate_sam_file} >> file_list.txt
        echo ~{flagstat_file} >> file_list.txt
        echo ~{star_log} >> file_list.txt

        qualimap_bamqc_dir=$(basename ~{qualimap_bamqc} ".tar.gz")
        tar -xzf ~{qualimap_bamqc}
        echo "$qualimap_bamqc_dir"/genome_results.txt >> file_list.txt
        for file in $(find "$qualimap_bamqc_dir"/raw_data_qualimapReport/); do
            echo $file >> file_list.txt
        done

        if [ "~{rnaseq}" = "true" ]; then
            qualimap_rnaseq_dir=$(basename ~{qualimap_rnaseq} ".tar.gz")
            tar -xzf ~{qualimap_rnaseq}
            echo "$qualimap_rnaseq_dir"/rnaseq_qc_results.txt >> file_list.txt
            for file in $(find "$qualimap_rnaseq_dir"/raw_data_qualimapReport/); do
                echo $file >> file_list.txt
            done
        else
            tar -xzf ~{fastq_screen_tar}
            fastq_screen_dir=$(basename ~{fastq_screen_tar} ".tar.gz")
            for file in $(find "$fastq_screen_dir"); do
                echo $file >> file_list.txt
            done
        fi

        tar -xzf ~{fastqc}
        fastqc_dir=$(basename ~{fastqc} ".tar.gz")
        for file in $(find "$fastqc_dir"); do
            echo $file >> file_list.txt
        done

        multiqc --cl_config "extra_fn_clean_exts: '_qualimap_bamqc_results'" \
            --file-list file_list.txt -o ~{out_directory}
        tar -czf ~{out_tar_gz} ~{out_directory}
    }

    runtime {
        disk: disk_size + " GB"
        docker: 'stjudecloud/multiqc:1.1.0'
        memory: memory_gb + " GB"
        maxRetries: max_retries
    }

    output {
        File out = "multiqc_results.tar.gz"
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "This WDL tool generates a MultiQC quality control metrics report summary from input QC result files."
    }

    parameter_meta {
        sorted_bam: "A aligned, sorted BAM file"
        validate_sam_file: "A file output from Picard's ValidateSam tool"
        qualimap_bamqc: "Tarballed directory of files output by Qualimap's BamQC mode"
        qualimap_rnaseq: "Tarballed directory of files output by Qualimap's RNA-seq mode"
        fastqc: "Tarballed directory of files output by FastQC"
        flagstat_file: "A file containing the output of Samtools' flagstat command for the input STAR aligned BAM file"
    }
}
