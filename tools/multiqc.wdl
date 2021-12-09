## # MultiQC
##
## This WDL tool wraps the [MultiQC](https://multiqc.info/) tool.
## MultiQC aggregates quality control results for bioinformatics.

version 1.0

task multiqc {
    input {
        File validate_sam_file
        File flagstat_file
        File qualimap_bamqc
        File? qualimap_rnaseq
        File fastqc
        File instrument_file
        File read_length_file
        File encoding_file
        File? strandedness_file
        File? junction_annotation
        File? fastq_screen
        File? star_log
        Int max_retries = 1
        Int memory_gb = 5
        Int disk_size = 20
    }

    String out_directory = basename(validate_sam_file, ".ValidateSamFile.txt") + ".multiqc"
    String out_tar_gz = out_directory + ".tar.gz"

    command {
        set -eo pipefail
        
        # set ENV variables for `multiqc`
        export LC_ALL=C.UTF-8
        export LANG=C.UTF-8
        
        (
            echo ~{validate_sam_file}
            echo ~{flagstat_file}
            echo ~{instrument_file}
            echo ~{read_length_file}
            echo ~{encoding_file}
            echo ~{star_log}
            echo ~{strandedness_file}
            echo ~{junction_annotation}
        ) > file_list.txt

        qualimap_bamqc_dir=$(basename ~{qualimap_bamqc} ".tar.gz")
        tar -xzf ~{qualimap_bamqc}
        echo "$qualimap_bamqc_dir"/genome_results.txt >> file_list.txt
        for file in "$qualimap_bamqc_dir"/raw_data_qualimapReport/*; do
            echo "$file" >> file_list.txt
        done

        tar -xzf ~{fastqc}
        fastqc_dir=$(basename ~{fastqc} ".tar.gz")
        for file in "$fastqc_dir"/*; do
            echo "$file" >> file_list.txt
        done
        
        if [ "~{if defined(qualimap_rnaseq) then "rnaseq" else ""}" = "rnaseq" ]; then
            qualimap_rnaseq_dir=$(basename ~{qualimap_rnaseq} ".tar.gz")
            tar -xzf ~{qualimap_rnaseq}
            echo "$qualimap_rnaseq_dir"/rnaseq_qc_results.txt >> file_list.txt
            echo "$qualimap_rnaseq_dir"/qualimap_counts.txt >> file_list.txt
            for file in "$qualimap_rnaseq_dir"/raw_data_qualimapReport/*; do
                echo "$file" >> file_list.txt
            done
        fi

        if [ "~{if defined(fastq_screen) then "wgs_or_wes" else ""}" = "wgs_or_wes" ]; then
            fastq_screen_dir=$(basename ~{fastq_screen} ".tar.gz")
            tar -xzf ~{fastq_screen}
            for file in "$fastq_screen_dir"/*; do
                echo "$file" >> file_list.txt
            done
        fi

        multiqc --verbose -c /home/.multiqc_config.yaml \
            --file-list file_list.txt -o ~{out_directory}
        
        tar -czf ~{out_tar_gz} ~{out_directory}
    }

    runtime {
        disk: disk_size + " GB"
        docker: 'ghcr.io/stjudecloud/multiqc:1.3.4'
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
        validate_sam_file: "A file output from Picard's ValidateSam tool"
        qualimap_bamqc: "Tarballed directory of files output by Qualimap's BamQC mode"
        qualimap_rnaseq: "Tarballed directory of files output by Qualimap's RNA-seq mode"
        fastqc: "Tarballed directory of files output by FastQC"
        fastq_screen: "Tarballed directory of files output by FastQ Screen"
        flagstat_file: "A file containing the output of Samtools' flagstat command for the input BAM file"
        star_log: "The log file of a STAR alignment run"
    }
}
