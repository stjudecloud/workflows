task print_version {
    command {
        fq --version
    }

    runtime {
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }

    output {
        String out = read_string(stdout())
    }

}

task fqlint {
    File read1
    File read2

    runtime {
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }

    command {
        fq lint ${read1} ${read2}
    }
}
