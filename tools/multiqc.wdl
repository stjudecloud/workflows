## Description: 
##
## This WDL tool wraps the MultiQC tool (https://multiqc.info/).
## MultiQC aggregates quality control results for bioinformatics.

task multiqc {
    File star
    String validate_sam_string
    Array[File] qualimap_bamqc
    Array[File] qualimap_rnaseq
    Array[File] fastqc_files
    File flagstat_file

    command {
        echo ${star} > file_list.txt
        echo ${validate_sam_string} > validate_sam.txt
        echo validate_sam.txt >> file_list.txt
        for file in ${sep=' ' qualimap_bamqc} ; do
            echo $file >> file_list.txt
        done
        for file in ${sep=' ' qualimap_rnaseq} ; do
            echo $file >> file_list.txt
        done
        for file in ${sep=' ' fastqc_files} ; do
            echo $file >> file_list.txt
        done
        echo ${flagstat_file} >> file_list.txt

        multiqc --file-list file_list.txt -o multiqc_results
    }

    runtime {
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }

    output {
        File out = "multiqc_results"
    }
    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
        author: "Andrew Frantz"
        email: "andrew.frantz@stjude.org"
        description: "This WDL tool generates a MultiQC quality control metrics report summary from input QC result files."
    }
    parameter_meta {
        star: "A STAR aligned BAM file"
        validate_sam_string: "A string output from Picard's ValidateSam tool"
        qualimap_bamqc: "An array of files output by Qualimap's BamQC mode"
        qualimap_rnaseq: "An array of files output by Qualimap's RNA-seq mode"
        fastqc_files: "An array of files output by FastQC"
        flagstat_file: "A file containing the output of Samtools' flagstat command for the input STAR aligned BAM file"
    }
}
