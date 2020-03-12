## Description:
##
## This WDL tool includes custom scripts to parse and validate QC output.  

version 1.0

task parse_validate_bam {
    input {
        String in
        Boolean strict = true
        Int max_retries = 1
    }
 
    command {
        if [ "~{strict}" == "true" ]
        then 
           if [ "$(echo "~{in}" | grep -c "ERROR")" -gt 0 ]
           then 
              echo "Errors detected by Picard ValidateSamFile" > /dev/stderr
              exit 1 
           fi
        elif [ "$(echo "~{in}" | grep -Ec "ERROR|WARNING")" -gt 0 ]
        then
              echo "Errors detected by Picard ValidateSamFile" > /dev/stderr
              exit 1 
        fi
    }

    runtime {
        docker: 'stjudecloud/util:1.0.0-alpha'
        maxRetries: max_retries
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "This WDL tool is a utility for parsing the output of Picard's ValidateSamFile command."
    }

    parameter_meta {
        in: "Output file from Picard ValidateSamFile"
    }
}