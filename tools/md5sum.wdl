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
        md5sum ${infile} > stdout.txt
    }
    output {
        String out = read_string("stdout.txt")
    }
}
task check_checksum {
    File infile
  
    command { 
        md5sum -c ${infile} > stdout.txt
    } 
    output {
        String out = read_string("stdout.txt")
    }
}
