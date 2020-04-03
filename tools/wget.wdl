## # wget 
##
## This WDL tool wraps the wget tool.

version 1.0

task wget_print_version {
    command {
        wget --version
    }

    runtime {
        docker: 'stjudecloud/util:1.0.0-alpha'
    }

    output {
        String out = read_string(stdout())
    }
}

task download {
    input {
        String url 
        String outfilename
        Int max_retries = 1
    }

    runtime {
        docker: 'stjudecloud/util:1.0.0-alpha'
        maxRetries: max_retries
    }

    command {
        wget ~{url} -O ~{outfilename}
    }

    output {
        File outfile = outfilename
    }

    meta {
        author: "Clay McLeod"
        email: "clay.mcleod@stjude.org"
        description: "This WDL tool uses wget to download a file from a remote URL to the local filesystem" 
    }

    parameter_meta {
        url: "URL of the file to download"
        outfilename: "Name to use for the output file"
    }
}
