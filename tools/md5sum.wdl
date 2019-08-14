task print_version {
    command {
        md5sum --version
    }

    runtime {
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
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

    runtime {
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
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

    runtime {
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }

    output {
        String out = read_string("stdout.txt")
    }
}
