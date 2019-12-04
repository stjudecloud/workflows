## Description:
##
## This WDL tool wraps the FastQC tool (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
## FastQC generates quality control metrics for sequencing pipelines. 

task fastqc {
    File bam
    Int ncpu
    String prefix = basename(bam, ".bam")

    Int bam_size = size(bam, "GiB")
    Int disk_size = ceil((bam_size * 2) + 10)

    command {
        mkdir ${prefix}_fastqc_results
        fastqc -f bam \
            -o ${prefix}_fastqc_results \
            -t ${ncpu} \
            ${bam}
    }

    runtime {
        disk: disk_size + " GB"
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }

    output {
        Array[File] out_files = glob("${prefix}_fastqc_results/*")
    }
}
