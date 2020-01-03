## Description: 
##
## This WDL tool wraps the SAMtools package (http://samtools.sourceforge.net/).
## SAMtools provides utlities for manipulating SAM format sequence alignments.

task print_version {
    command {
        samtools --version
    }

    runtime {
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }

    output {
        String out = read_string(stdout())
    }

}

task quickcheck {
    File bam

    command {
        samtools quickcheck ${bam}
    }

    runtime {
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }
    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
        author: "Andrew Frantz"
        email: "andrew.frantz@stjude.org"
        description: "This WDL tool runs Samtools quickcheck on the input BAM file. This checks that the BAM file appears to be intact, e.g. header exists, at least one sequence is present, and the end-of-file marker exists."
    }
    parameter_meta {
        bam: "Input BAM format file to generate coverage for"
    }
}

task split {
    File bam
    Boolean? reject_unaccounted
    String prefix = basename(bam, ".bam")

    command {
        samtools split -u ${prefix}.unaccounted_reads.bam -f '%*_%!.%.' ${bam}
        samtools view ${prefix}.unaccounted_reads.bam > unaccounted_reads.bam
        if [ ${default='true' reject_unaccounted} -a -s unaccounted_reads.bam ] 
            then exit 1; 
            else rm ${prefix}.unaccounted_reads.bam
        fi 
        rm unaccounted_reads.bam
    }
 
    runtime {
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }
   
    output {
       Array[File] split_bams = glob("*.bam")
    }
    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
        author: "Andrew Frantz"
        email: "andrew.frantz@stjude.org"
        description: "This WDL tool runs Samtools split on the input BAM file. This splits the BAM by read group into one or more output files. It optionally errors if there are reads present that do not belong to a read group."
    }
    parameter_meta {
        bam: "Input BAM format file to generate coverage for"
        reject_unaccounted: "If true, error if there are reads present that do not have read group information."
    }
}

task flagstat { 
    File bam

    String outfile = basename(bam, ".bam")+".flagstat.txt"

    command {
        samtools flagstat ${bam} > ${outfile}
    }

    runtime {
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }

    output { 
       File flagstat = outfile
    }
    meta {
        authors: [
           {
              name: "Andrew Thrasher"
              email: "andrew.thrasher@stjude.org"
           }, {
              name: "Andrew Frantz"
              email: "andrew.frantz@stjude.org"
           }
        ]
        description: "This WDL tool generates a FastQC quality control metrics report for the input BAM file."
    }
    parameter_meta {
        bam: "Input BAM format file to generate coverage for"
    }
}

task index {
    File bam
    
    String name = basename(bam)
    String outfile = basename(bam)+".bai"

    command {
        samtools index ${bam} ${outfile}
    }

    runtime {
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }

    output {
       File bai = outfile
    }
    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
        author: "Andrew Frantz"
        email: "andrew.frantz@stjude.org"
        description: "This WDL tool runs Samtools flagstat on the input BAM file. Produces statistics about the alignments based on the bit flags set in the BAM."
    }
    parameter_meta {
        bam: "Input BAM format file to generate coverage for"
    }
}
