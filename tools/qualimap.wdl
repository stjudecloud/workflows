## Description: 
##
## This WDL tool wraps the QualiMap tool (http://qualimap.bioinfo.cipf.es/).
## QualiMap computes metrics to facilitate evaluation of sequencing data. 

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
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }

    output {
        Array[File] out_files = glob("${prefix}_qualimap_results/*")
    }
    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
        author: "Andrew Frantz"
        email: "andrew.frantz@stjude.org"
        description: "This WDL tool runs QualiMap's bamqc tool on the input BAM file."
    }
    parameter_meta {
        bam: "Input BAM format file to generate coverage for"
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
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }

    output {
        Array[File] out_files = glob("${outdir}/*")
    }
    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
        author: "Andrew Frantz"
        email: "andrew.frantz@stjude.org"
        description: "This WDL tool generates runs QualiMap's rnaseq tool on the input BAM file."
    }
    parameter_meta {
        bam: "Input BAM format file to generate coverage for"
        gencode_gtf: "A GTF format features file containing Gencode features"
        strand: "Strand information for RNA-seq experiments. Options: [strand-specific-forward, strand-specific-reverse, non-strand-specific]"
    }
}
