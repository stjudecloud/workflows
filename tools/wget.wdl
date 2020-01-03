## Description: 
##
## This WDL tool wraps the wget tool.

task print_version {
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
    String url 
    String outfilename

    runtime {
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }

    command {
      wget ${url} -O ${outfilename}
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
