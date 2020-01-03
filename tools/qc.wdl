## Description: 
##
## This WDL tool includes custom scripts to parse and validate QC output.  

task parse_validate_bam {
    String in
    Boolean strict = true
 
    command {
        if [ "${strict}" == "true" ]
        then 
           if [ $(echo "${in}" | grep -c "ERROR") -gt 0 ] 
           then 
              echo "Errors detected by Picard ValidateSamFile" > /dev/stderr
              exit -1 
           fi
        elif [ $(echo "${in}" | grep -Ec "ERROR|WARNING") -gt 0 ]
        then
              echo "Errors detected by Picard ValidateSamFile" > /dev/stderr
              exit -1 
        fi
    }

    runtime {
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
    }
    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
        author: "Andrew Frantz"
        email: "andrew.frantz@stjude.org"
        description: "This WDL tool is a utility for parsing the output of Picard's ValidateSamFile command."
    }
    parameter_meta {
        bam: "Input BAM format file to generate coverage for"
    }
}

task parse_infer_experiment {
    String in

    command { 
        if [ $(echo "${in}" | grep -c 'PairEnd') -gt 0 ]
        then
            if (( $(echo "$(echo "${in}" | tail -n 1 | sed  's/^.*: //') > 0.5" | bc -l ) ))
            then
               echo 'strand-specific-reverse' > stdout.txt
            elif (( $(echo "$(echo "${in}" | tail -n 2 | head -n 1 | sed  's/^.*: //') > 0.5" | bc -l ) ))
            then 
               echo 'strand-specific-forward' > stdout.txt
            else
               echo 'non-strand-specific' > stdout.txt
            fi
        elif [ $(echo "${in}" | grep -c 'SingleEnd') -gt 0 ]
        then 
            if (( $(echo "$(echo "${in}" | tail -n 1 | sed  's/^.*: //') > 0.5" | bc -l ) ))
            then
               echo 'strand-specific-reverse' > stdout.txt
            elif (( $(echo "$(echo "${in}" | tail -n 2 | head -n 1 | sed  's/^.*: //') > 0.5" | bc -l ) ))
            then 
               echo 'strand-specific-forward' > stdout.txt
            else
               echo 'non-strand-specific' > stdout.txt
            fi
        else
           echo "infer_experiment failed to determine type" > /dev/stderr
           exit -1 
        fi
    }

    runtime {
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
        description: "This WDL tool parses the output of RSeQC's infer_experiment package."
    }
    parameter_meta {
        bam: "Input BAM format file to generate coverage for"
    }
}
