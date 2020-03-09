## Description:
##
## This WDL tool wraps the MultiQC tool (https://multiqc.info/).
## MultiQC aggregates quality control results for bioinformatics.

version 1.0

task multiqc {
    input {
        File sorted_bam
        String validate_sam_string
        File qualimap_bamqc
        File? qualimap_rnaseq
        Array[File] fastqc_files
        Array[File] fastq_screen
        File flagstat_file
        File? bigwig_file
        File? star_log
        Int max_retries = 1
        String wait_var = ""
    }

    Float star_size = size(sorted_bam, "GiB")
    Int disk_size = ceil((star_size * 4) + 10)

    command {
        echo ~{sorted_bam} > file_list.txt
        echo ~{validate_sam_string} > validate_sam.txt
        echo validate_sam.txt >> file_list.txt

        mkdir qualimap_bamqc/
        tar -xzf ~{qualimap_bamqc} --strip-components=1 -C qualimap_bamqc/;
        echo qualimap_bamqc/genome_results.txt >> file_list.txt
        for file in `find qualimap_bamqc/raw_data_qualimapReport/`; do
            echo $file >> file_list.txt
        done

        mkdir qualimap_rnaseq/
        tar -xzf ~{qualimap_rnaseq} --strip-components=1 -C qualimap_rnaseq/;
        echo qualimap_rnaseq/rnaseq_qc_results.txt >> file_list.txt
        for file in `find qualimap_rnaseq/raw_data_qualimapReport/`; do
            echo $file >> file_list.txt
        done

        for file in ~{sep=' ' fastqc_files} ; do
            echo $file >> file_list.txt
        done

        for file in ~{sep=' ' fastq_screen} ; do
            echo $file >> file_list.txt
        done

        # shellcheck disable=SC2129
        echo ~{flagstat_file} >> file_list.txt
        echo ~{bigwig_file} >> file_list.txt
        echo ~{star_log} >> file_list.txt

        multiqc --file-list file_list.txt -o multiqc_results
        tar -czf multiqc_results.tar.gz multiqc_results
    }

    runtime {
        disk: disk_size + " GB"
        docker: 'stjudecloud/multiqc:1.0.0-alpha'
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
        validate_sam_string: "A string output from Picard's ValidateSam tool"
        qualimap_bamqc: "An array of files output by Qualimap's BamQC mode"
        qualimap_rnaseq: "An array of files output by Qualimap's RNA-seq mode"
        fastqc_files: "An array of files output by FastQC"
        flagstat_file: "A file containing the output of Samtools' flagstat command for the input STAR aligned BAM file"
    }
}
