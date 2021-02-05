## # MultiQC
##
## This WDL tool wraps the [MultiQC](https://multiqc.info/) tool.
## MultiQC aggregates quality control results for bioinformatics.

version 1.0

task multiqc {
    input {
        File bam
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

    String out_directory = basename(bam, ".bam") + "_multiqc"
    String out_tar_gz = out_directory + ".tar.gz"
    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil((bam_size * 2) + 10)

    command {
        set -eo pipefail
        
        # set ENV variables for `multiqc`
        export LC_ALL=C.UTF-8
        export LANG=C.UTF-8
        
        echo ~{bam} > file_list.txt
        echo ~{validate_sam_file} >> file_list.txt
        echo ~{flagstat_file} >> file_list.txt

        qualimap_bamqc_dir=$(basename ~{qualimap_bamqc} ".tar.gz")
        tar -xzf ~{qualimap_bamqc}
        echo "$qualimap_bamqc_dir"/genome_results.txt >> file_list.txt
        for file in $(find "$qualimap_bamqc_dir"/raw_data_qualimapReport/); do
            echo "$file" >> file_list.txt
        done

        if [ "~{if defined(qualimap_rnaseq) then "rnaseq" else ""}" = "rnaseq" ]; then
            echo ~{star_log} >> file_list.txt
            qualimap_rnaseq_dir=$(basename ~{qualimap_rnaseq} ".tar.gz")
            tar -xzf ~{qualimap_rnaseq}
            echo "$qualimap_rnaseq_dir"/rnaseq_qc_results.txt >> file_list.txt
            for file in $(find "$qualimap_rnaseq_dir"/raw_data_qualimapReport/); do
                echo "$file" >> file_list.txt
            done
        fi
        if [ "~{if defined(fastq_screen) then "wgs_or_wes" else ""}" = "wgs_or_wes" ]; then
            tar -xzf ~{fastq_screen}
            fastq_screen_dir=$(basename ~{fastq_screen} ".tar.gz")
            for file in $(find "$fastq_screen_dir"); do
                echo "$file" >> file_list.txt
            done
        fi

        tar -xzf ~{fastqc}
        fastqc_dir=$(basename ~{fastqc} ".tar.gz")
        for file in $(find "$fastqc_dir"); do
            echo "$file" >> file_list.txt
        done

        multiqc --cl_config "extra_fn_clean_exts: '_qualimap_bamqc_results'" \
            --file-list file_list.txt -o ~{out_directory}
        
        tar -czf ~{out_tar_gz} ~{out_directory}
    }

    runtime {
        disk: disk_size + " GB"
        docker: 'stjudecloud/multiqc:1.2.0'
        memory: memory_gb + " GB"
        maxRetries: max_retries
    }

    output {
        File out = out_tar_gz
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "This WDL tool generates a MultiQC quality control metrics report summary from input QC result files."
    }

    parameter_meta {
        bam: "The BAM referenced by other input reports"
        validate_sam_file: "A file output from Picard's ValidateSam tool"
        qualimap_bamqc: "Tarballed directory of files output by Qualimap's BamQC mode"
        qualimap_rnaseq: "Tarballed directory of files output by Qualimap's RNA-seq mode"
        fastqc: "Tarballed directory of files output by FastQC"
        fastq_screen: "Tarballed directory of files output by FastQ Screen"
        flagstat_file: "A file containing the output of Samtools' flagstat command for the input BAM file"
        star_log: "The log file of a STAR alignment run"
    }
}
