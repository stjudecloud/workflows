task multiqc {
    File star
    File dups
    String validate_sam_string
    Array[File] qualimap_bamqc
    Array[File] qualimap_rnaseq
    Array[File] fastqc
    File flagstat

    command {
        echo ${star} > file_list.txt
        echo ${dups} >> file_list.txt
        echo ${validate_sam_string} > validate_sam.txt
        echo validate_sam.txt >> file_list.txt
        for file in ${sep=' ' qualimap_bamqc} ; do
            echo $file >> file_list.txt
        done
        for file in ${sep=' ' qualimap_rnaseq} ; do
            echo $file >> file_list.txt
        done
        for file in ${sep=' ' fastqc} ; do
            echo $file >> file_list.txt
        done
        echo ${flagstat} >> file_list.txt

        multiqc --file-list file_list.txt -o multiqc_results
    }

    output {
        File out = "multiqc_results"
    }
}