## Description:
##
## This WDL tool includes custom scripts to parse and validate QC output.  

version 1.0

task parse_infer_experiment {
    input {
        String in
        Int max_retries = 1
    }

    command { 
        if [ "$(echo "~{in}" | grep -c 'PairEnd')" -gt 0 ]
        then
            if (( $(echo "$(echo "~{in}" | tail -n 1 | sed  's/^.*: //') > 0.5" | bc -l ) ))
            then
               echo 'strand-specific-reverse' > stdout.txt
            elif (( $(echo "$(echo "~{in}" | tail -n 2 | head -n 1 | sed  's/^.*: //') > 0.5" | bc -l ) ))
            then 
               echo 'strand-specific-forward' > stdout.txt
            else
               echo 'non-strand-specific' > stdout.txt
            fi
        elif [ "$(echo "~{in}" | grep -c 'SingleEnd')" -gt 0 ]
        then 
            if (( $(echo "$(echo "~{in}" | tail -n 1 | sed  's/^.*: //') > 0.5" | bc -l ) ))
            then
               echo 'strand-specific-reverse' > stdout.txt
            elif (( $(echo "$(echo "~{in}" | tail -n 2 | head -n 1 | sed  's/^.*: //') > 0.5" | bc -l ) ))
            then 
               echo 'strand-specific-forward' > stdout.txt
            else
               echo 'non-strand-specific' > stdout.txt
            fi
        else
           echo "infer_experiment failed to determine type" > /dev/stderr
           exit 1 
        fi
    }

    runtime {
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
        maxRetries: max_retries
    }

    output {
        String out = read_string("stdout.txt")
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "This WDL tool parses the output of RSeQC's infer_experiment.py script."
    }

    parameter_meta {
        in: "Output file from RSeQC's infer_experiment.py script."
    }
}
