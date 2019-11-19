## Description: 
##
## This WDL tool wraps any functionality related to gzip.

task print_version {
    command {
        gzip --version
    }

    runtime {
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }

    output {
        String out = read_string(stdout())
    }

}

task unzip {
    File infile
    String outfilename = basename(infile, ".gz")

    runtime {
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }

    command {
			gunzip ${infile} -c > ${outfilename}
    }

		output {
			File outfile = outfilename
		}
}