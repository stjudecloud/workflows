task print_version {
    command {
        md5sum --version
    }

    output {
        String out = read_string(stdout())
    }

}
