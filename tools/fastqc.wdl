## Description: 
##
## This WDL tool wraps the FastQC tool (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
## FastQC generates quality control metrics for sequencing pipelines. 

task fastqc {
    File bam
    Int ncpu
    String prefix = basename(bam, ".bam")

    command {
        mkdir ${prefix}_fastqc_results
        fastqc -f bam \
            -o ${prefix}_fastqc_results \
            -t ${ncpu} \
            ${bam}
    }

    runtime {
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }

    output {
        Array[File] out_files = glob("${prefix}_fastqc_results/*")
    }
    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
        author: "Andrew Frantz"
        email: "andrew.frantz@stjude.org"
        description: "This WDL tool generates a FastQC quality control metrics report for the input BAM file."
    }
    parameter_meta {
        bam: "Input BAM format file to generate coverage for"
    }
}
