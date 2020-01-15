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

    runtime {
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }

    command {
        fq lint ${read1} ${read2}
    }
    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
        author: "Andrew Frantz"
        email: "andrew.frantz@stjude.org"
        description: "This WDL tool performs quality control on the input FastQ pairs to ensure proper formatting."
    }
    parameter_meta {
        read1: "Input FastQ with read one"
        read2: "Input FastQ with read two"
    }
}
