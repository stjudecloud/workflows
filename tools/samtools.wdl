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
}
