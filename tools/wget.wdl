## Description: 
##
## This WDL tool wraps the wget tool.

version 1.0

task wget_print_version {
    command {
        wget --version
    }

    runtime {
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }

    output {
        String out = read_string(stdout())
    }

}

task download {
    input {
        String url 
        String outfilename
    }

    runtime {
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }

    command {
      wget ${url} -O ${outfilename}
    }

    output {
      File outfile = outfilename
    }
}
