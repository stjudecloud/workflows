task print_version {
    command {
        samtools --version
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
}

task split {
    File bam
    Boolean? reject_unaccounted
    String basename = basename(bam, ".bam")

    command {
        samtools split -u ${basename}.unaccounted_reads.bam -f '%*_%!.%.' ${bam}; 
        samtools view ${basename}.unaccounted_reads.bam > unaccounted_reads.bam; 
        if [ ${default='true' reject_unaccounted} -a -s unaccounted_reads.bam ] 
            then exit 1; 
            else rm ${basename}.unaccounted_reads.bam
        fi 
        rm unaccounted_reads.bam
    }
    
    output {
       Array[File] out_bams = glob("*.bam")
    }
}

task flagstat { 
    File bam

    String outfile = basename(bam, ".bam")+".flagstat.txt"

    command {
        samtools flagstat ${bam} > ${outfile}
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
    output {
       File bai = outfile
    }
}
