task bamqc {
    File bam
    Int ncpu
    String basename = basename(bam, ".bam")

    command {
        qualimap bamqc -bam ${bam} \
            -outdir ${basename}_qualimap_results \
            -nt ${ncpu}; \
        cd ${basename}_fastqc_results; \
        dir=`pwd`; \
        ln *.gz "$dir"; \
        ls > file-list.txt
    }

    output {
        Array[File] files = read_lines("file-list.txt")
    }
}
