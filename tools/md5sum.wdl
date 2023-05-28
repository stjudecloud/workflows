## Description:
##
## This WDL file wraps the md5sum tool from the [GNU core utilities](https://github.com/coreutils/coreutils).
## md5sum is a utility for generating and verifying MD5 hashes.

version 1.0

task compute_checksum {
    meta {
        description: "This WDL task generates an MD5 checksum for the input file."
    }

    parameter_meta {
        infile: "Input file to generate MD5 checksum for"
    }

    input {
        File infile
        String outfile_name = basename(infile) + ".md5"
        Int memory_gb = 5
        Int max_retries = 1
    }

    Float infile_size = size(infile, "GiB")
    Int disk_size_gb = ceil(infile_size * 1.5)

    command <<<
        md5sum ~{infile} > ~{outfile_name}
    >>>

    output {
        File md5sum = outfile_name
    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        docker: 'ghcr.io/stjudecloud/util:1.2.0'
        maxRetries: max_retries
    }
}
