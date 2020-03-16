## Description:
##
## This WDL tool includes custom scripts to parse and reformat 
## task output as part of a workflow. 

version 1.0

task get_read_groups {
    input {
        File bam
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil((bam_size * 2) + 10)

    command <<<
        samtools view -H ~{bam} | grep "@RG" \
            | cut -f 2- \
            | sed -e 's/\t/ /g' \
            | awk '{print}' ORS=' , ' \
            | sed 's/ , $//' >> stdout.txt
    >>>

    runtime {
        disk: disk_size + " GB"
        docker: 'stjudecloud/bioinformatics-base:bleeding-edge'
        maxRetries: max_retries
    }

    output { 
        String out = read_string("stdout.txt")
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "This WDL tool is a utility to get read group information from a BAM file and write it out to as a string" 
    }

    parameter_meta {
        bam: "Input BAM format file to generate coverage for"
    }
}

task file_prefix {
    input {
        File in_file
        Int max_retries = 1
    }

    command <<<
        basename ~{in_file} | awk -F '.' '{print $1}' > stdout.txt
    >>>

    runtime {
        docker: 'stjudecloud/util:1.0.0-alpha'
        maxRetries: max_retries
    }

    output { 
        String out = read_string("stdout.txt")
    }
}
