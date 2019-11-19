## Description: 
##
## This WDL tool wraps the wget tool.

task print_version {
    command {
        wget --version
    }

    runtime {
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }

    output {
        String out = read_string(stdout())
    }

}

task download {
    String url 
    String outfilename

    runtime {
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }

    command {
			wget ${url} -O ${outfilename}
    }

		output {
			File outfile = outfilename
		}
}
