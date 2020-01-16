## Description:
##
## This WDL tool wraps the md5sum tool from the GNU core 
## utilities (https://github.com/coreutils/coreutils).
## md5sum is a utility for generating and verifying MD5
## hashes.  

version 1.0

task md5sum_print_version {
    command {
        md5sum --version
    }

    runtime {
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }

    output {
        String out = read_string(stdout())
    }
}
task compute_checksum {
    input {
        File infile
    }

    String outfile = infile + ".md5"
    Float infile_size = size(infile, "GiB")
    Int disk_size = ceil((infile_size * 2) + 10)

    command {
        md5sum ${infile} > ${outfile}
    }

    runtime {
        disk: disk_size + " GB"
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }

    output {
        String stdout = read_string(outfile)
        File outfile = outfile
    }
    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
        author: "Andrew Frantz"
        email: "andrew.frantz@stjude.org"
        description: "This WDL tool generates an MD5 checksum for the input file."
    }
    parameter_meta {
        infile: "Input file to generate MD5 checksum"
    }
}
task check_checksum {
    input {
        File infile
    }
    Float infile_size = size(infile, "GiB")
    Int disk_size = ceil((infile_size * 2) + 10)

    command { 
        md5sum -c ${infile} > stdout.txt
    } 

    runtime {
        disk: disk_size + " GB"
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }

    output {
        String out = read_string("stdout.txt")
    }
    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
        author: "Andrew Frantz"
        email: "andrew.frantz@stjude.org"
        description: "This WDL tool checks a list of MD5 checksums against the corresponding files to verify integrity" 
    }
    parameter_meta {
        infile: "Input file containing checksums to check" 
    }
}
