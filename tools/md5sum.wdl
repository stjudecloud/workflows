task print_version {
    command {
        md5sum --version
    }

    output {
        String out = read_string(stdout())
    }
}
task compute_checksum {
    File infile 

    command {
        md5sum ${infile}
    }
    output {
        String out = read_string(stdout())
    }
}
task check_checksum {
    File infile
  
    command { 
        md5sum -c ${infile}
    } 
    output {
        String out = read_string(stdout())
    }
}
