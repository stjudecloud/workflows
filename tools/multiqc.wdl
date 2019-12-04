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

    Int star_size = size(star, "GiB")
    Int disk_size = ceil((bam_size * 4) + 10)

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
        zip -r multiqc_results.zip multiqc_results
    }

    runtime {
        disk: disk_size + " GB"
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }

    output {
        File out = "multiqc_results.zip"
    }
}
