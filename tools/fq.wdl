## Description:
##
## This WDL tool wraps the fq tool (https://github.com/stjude/fqlib).
## The fq library provides methods for manipulating Illumina generated 
## FastQ files.

task print_version {
    command {
        fq --version
    }

    runtime {
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }

    output {
        String out = read_string(stdout())
    }
}

task fqlint {
    File read1
    File read2

    Float read1_size = size(read1, "GiB")
    Float read2_size = size(read2, "GiB")
    Int disk_size = ceil(((read1_size + read2_size) * 2) + 10)

    runtime {
        disk: disk_size + " GB"
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }

    command {
        fq lint ${read1} ${read2}
    }
}
