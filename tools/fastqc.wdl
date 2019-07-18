task fastqc {
    File bam
    Int ncpu
    String basename = basename(bam, ".bam")

    command {
        fastqc -f bam \
            -o ${basename}_fastqc_results \
            -t ${ncpu} \
            ${bam};
    }

    output {
        Array[File] out_files = glob("${basename}_fastqc_results/*")
    }
}
