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
        file: "Input file to generate MD5 checksum for"  # TODO rename
    }

    input {
        File file
        String outfile_name = basename(file) + ".md5"
        Int memory_gb = 5
        Int modify_disk_size_gb = 0
        Int max_retries = 1
    }

    Float file_size = size(file, "GiB")
    Int disk_size_gb = ceil(file_size) + 10 + modify_disk_size_gb

    command <<<
        md5sum ~{file} > ~{outfile_name}
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
