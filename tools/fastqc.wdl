task fastqc {
    File bam
    Int ncpu
    String basename = basename(bam, ".bam")

    command {
        mkdir ${basename}_fastqc_results
        fastqc -f bam \
            -o ${basename}_fastqc_results \
            -t ${ncpu} \
            ${bam}
    }

    output {
        Array[File] out_files = glob("${basename}_fastqc_results/*")
    }
}
