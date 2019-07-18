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
    String? unaccounted_bam

    command {
        samtools split -u ${default="unaccounted.bam" unaccounted_bam} -f '%*_%!.%.' ${bam}
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
