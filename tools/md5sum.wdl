## [Homepage](https://github.com/coreutils/coreutils)

version 1.1

task compute_checksum {
    meta {
        description: "Generates an MD5 checksum for the input file"
        outputs: {
            md5sum: "STDOUT of the `md5sum` command that has been redirected to a file"
        }
    }

    parameter_meta {
        file: "Input file to generate MD5 checksum for"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File file
        Int modify_disk_size_gb = 0
    }

    Float file_size = size(file, "GiB")
    Int disk_size_gb = ceil(file_size) + 10 + modify_disk_size_gb

    String outfile_name = basename(file) + ".md5"

    command <<<
        md5sum "~{file}" > "~{outfile_name}"
    >>>

    output {
        File md5sum = outfile_name
    }

    runtime {
        memory: "4 GB"
        disks: "~{disk_size_gb} GB"
        container: "ghcr.io/stjudecloud/util:branch-shellcheck-2.2.1"
        maxRetries: 1
    }
}
