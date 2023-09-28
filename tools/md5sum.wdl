## Description:
##
## This WDL file wraps the md5sum tool from the [GNU core utilities](https://github.com/coreutils/coreutils).
## md5sum is a utility for generating and verifying MD5 hashes.

version 1.1

task compute_checksum {
    meta {
        description: "This WDL task generates an MD5 checksum for the input file."
        outputs: {
            md5sum: "STDOUT of the `md5sum` command that has been redirected to a file"
        }
    }

    parameter_meta {
        file: "Input file to generate MD5 checksum for"
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
        max_retries: "Number of times to retry in case of failure"
    }

    input {
        File file
        Int memory_gb = 4
        Int modify_disk_size_gb = 0
        Int max_retries = 1
    }

    Float file_size = size(file, "GiB")
    Int disk_size_gb = ceil(file_size) + 10 + modify_disk_size_gb

    String outfile_name = basename(file) + ".md5"

    command <<<
        md5sum ~{file} > ~{outfile_name}
    >>>

    output {
        File md5sum = outfile_name
    }

    runtime {
        memory: "~{memory_gb} GB"
        disk: "~{disk_size_gb} GB"
        docker: 'ghcr.io/stjudecloud/util:1.3.0'
        maxRetries: max_retries
    }
}
