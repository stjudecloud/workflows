task print_version {
    command {
        fq --version
    }

    output {
        String out = read_string(stdout())
    }

}
