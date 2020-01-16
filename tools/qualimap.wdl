## Description:
##
## This WDL tool wraps the QualiMap tool (http://qualimap.bioinfo.cipf.es/).
## QualiMap computes metrics to facilitate evaluation of sequencing data. 

version 1.0

task bamqc {
    input {
        File bam
        Int ncpu = 1
        String prefix = basename(bam, ".bam")
    }
    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil((bam_size * 2) + 10)

    command {
        qualimap bamqc -bam ${bam} \
            -outdir ${prefix}_qualimap_results \
            -nt ${ncpu} \
            --java-mem-size=6g
    }

    runtime {
        memory: "8 GB"
        disk: disk_size + " GB"
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }

    output {
        Array[File] out_files = glob("${prefix}_qualimap_results/*")
    }
}

task rnaseq {
    input {
        File bam
        File gencode_gtf
        String outdir = "qualimap_rnaseq"
        String strand = "strand-specific-reverse"
    }
    Float bam_size = size(bam, "GiB")
    Float gencode_gtf_size = size(gencode_gtf, "GiB")
    Int disk_size = ceil(((bam_size + gencode_gtf_size) * 6) + 10)
 
    command {
        qualimap rnaseq -bam ${bam} -gtf ${gencode_gtf} -outdir ${outdir} -oc qualimap_counts.txt -p ${strand} -pe --java-mem-size=14G
    }

    runtime {
        memory: "16 GB"
        disk: disk_size + " GB"
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }

    output {
        Array[File] out_files = glob("${outdir}/*")
    }
}
