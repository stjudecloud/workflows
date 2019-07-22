task bamqc {
    File bam
    Int ncpu
    String basename = basename(bam, ".bam")

    command {
        qualimap bamqc -bam ${bam} \
            -outdir ${basename}_qualimap_results \
            -nt ${ncpu}
    }

    output {
        Array[File] out_files = glob("${basename}_qualimap_results/*")
    }
}
