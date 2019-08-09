task bamqc {
    File bam
    Int ncpu
    String prefix = basename(bam, ".bam")

    command {
        qualimap bamqc -bam ${bam} \
            -outdir ${prefix}_qualimap_results \
            -nt ${ncpu} \
            --java-mem-size=6g
    }

    runtime {
        memory: "8G"
    }

    output {
        Array[File] out_files = glob("${prefix}_qualimap_results/*")
    }
}
task rnaseq {
    File bam
    File gencode_gtf
    String outdir = "qualimap_rnaseq"
    String strand = "strand-specific-reverse"
 
    command {
        qualimap rnaseq -bam ${bam} -gtf ${gencode_gtf} -outdir ${outdir} -oc qualimap_counts.txt -p ${strand} -pe --java-mem-size=6G
    }

    runtime {
        memory: "8G"
    }

    output {
        Array[File] out_files = glob("${outdir}/*")
    }
}
