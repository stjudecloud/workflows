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

}

task flagstat { 

}

task index {

}

