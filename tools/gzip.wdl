## Description: 
##
## This WDL tool wraps any functionality related to gzip.

version 1.0

task gzip_print_version {
    command {
        gzip --version
    }

    runtime {
        docker: 'docker.pkg.github.com/stjudecloud/workflows/util:0.1.0'
    }

    output {
        String out = read_string(stdout())
    }
}

task unzip {
    input {
        File infile
        String outfilename = basename(infile, ".gz")
        Int max_retries = 1
    }

    runtime {
        docker: 'docker.pkg.github.com/stjudecloud/workflows/util:0.1.0'
        maxRetries: max_retries
    }

    command {
        gunzip ${infile} -c > ${outfilename}
    }

    output {
        File outfile = outfilename
    }

    meta {
        author: "Clay McLeod"
        email: "clay.mcleod@stjude.org"
        description: "This WDL tool unzips a gzip archive" 
    }

    parameter_meta {
        infile: "Archive in gzip format to be extracted"
    }
}
