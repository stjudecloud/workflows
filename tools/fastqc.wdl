task fastqc {
    File bam
    Int ncpu
    String basename = basename(bam, ".bam")

    command {
        fastqc -f bam \
            -o ${basename}_fastqc_results \
            -t ${ncpu} \
            ${bam}; \
        cd ${basename}_fastqc_results; \
        dir=`pwd`; \
        ln *.gz "$dir"; \
        ls > file-list.txt
    }

    output {
        Array[File] files = read_lines("file-list.txt")
    }
}
