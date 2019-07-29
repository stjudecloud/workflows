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
task rnaseq {
    File bam
    File gencode_gtf
    String outdir = "qualimap_rnaseq"
    String strand = "strand-specific-reverse"
 
    command {
        qualimap rnaseq -bam ${bam} -gtf ${gencode_gtf} -outdir ${outdir} -oc qualimap_counts.txt -p ${strand} -pe
    }

    output {
        Array[File] out_files = glob("${outdir}/*")
    }
}
