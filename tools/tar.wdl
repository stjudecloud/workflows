## # tar 
##
## This WDL tool wraps any functionality related to tar.

version 1.0

task tar_print_version {
    command {
        tar --version
    }

    runtime {
        docker: 'stjudecloud/util:1.0.0'
    }

    output {
        String out = read_string(stdout())
    }
}

task untar {
    input {
        File infile
        Int max_retries = 1
    }

    runtime {
        docker: 'stjudecloud/util:1.0.0'
        maxRetries: max_retries
    }

    command {
        mkdir tar_output
        tar --no-same-owner -xf ${infile} -C tar_output
    }

    output {
        Array[File] outfiles = glob("tar_output/*")
    }

    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
        description: "This WDL tool unzips a tar archive" 
    }

    parameter_meta {
        infile: "Archive in tar format to be extracted"
    }
}
