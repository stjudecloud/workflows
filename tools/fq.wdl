task print_version {
    command {
        fq --version
    }

    output {
        String out = read_string(stdout())
    }

}

task fqlint {
    File read1
    File read2

    command {
        fq lint ${read1} ${read2}
    }
}