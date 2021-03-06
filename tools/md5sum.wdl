## Description:
##
## This WDL tool wraps the md5sum tool from the [GNU core utilities](https://github.com/coreutils/coreutils).
## md5sum is a utility for generating and verifying MD5
## hashes.  

version 1.0

task compute_checksum {
    input {
        File infile
        Int max_retries = 1
        Int memory_gb = 5
    }

    String outfilename = basename(infile) + ".md5"
    Float infile_size = size(infile, "GiB")
    Int disk_size = ceil((infile_size * 2) + 10)

    command {
        md5sum ~{infile} > ~{outfilename}
    }

    runtime {
        disk: disk_size + " GB"
        memory: memory_gb + " GB"
        docker: 'stjudecloud/util:1.0.0'
        maxRetries: max_retries
    }

    output {
        File outfile = outfilename
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "This WDL tool generates an MD5 checksum for the input file."
    }

    parameter_meta {
        infile: "Input file to generate MD5 checksum"
    }
}

task check_checksum {
    input {
        File infile
        Int max_retries = 1
    }

    String outfilename = basename(infile) + ".md5_check"
    Float infile_size = size(infile, "GiB")
    Int disk_size = ceil((infile_size * 2) + 10)

    command { 
        md5sum -c ~{infile} > ~{outfilename}
    } 

    runtime {
        disk: disk_size + " GB"
        docker: 'stjudecloud/util:1.0.0'
        maxRetries: max_retries
    }

    output {
        File out = outfilename
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "This WDL tool checks a list of MD5 checksums against the corresponding files to verify integrity" 
    }

    parameter_meta {
        infile: "Input file containing checksums to check" 
    }
}
