## Description:
##
## This WDL file wraps the md5sum tool from the [GNU core utilities](https://github.com/coreutils/coreutils).
## md5sum is a utility for generating and verifying MD5
## hashes.  

version 1.0

task compute_checksum {
    input {
        File infile
        String outfile_name = basename(infile) + ".md5"
        Int memory_gb = 5
        Int max_retries = 1
    }

    Float infile_size = size(infile, "GiB")
    Int disk_size = ceil((infile_size * 2) + 10)

    command {
        md5sum ~{infile} > ~{outfile_name}
    }

    runtime {
        disk: disk_size + " GB"
        memory: memory_gb + " GB"
        docker: 'ghcr.io/stjudecloud/util:1.2.0'
        maxRetries: max_retries
    }

    output {
        File md5sum = outfile_name
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "This WDL task generates an MD5 checksum for the input file."
    }

    parameter_meta {
        infile: "Input file to generate MD5 checksum"
    }
}
