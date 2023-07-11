## # MultiQC
##
## This WDL file wraps the [MultiQC](https://multiqc.info/) tool.
## MultiQC aggregates quality control results for bioinformatics.

version 1.1

task multiqc {
    meta {
        description: "This WDL task generates a MultiQC quality control metrics report summary from input QC result files."
    }

    parameter_meta {
        input_files: "An array of files for MultiQC to compile into a report. Invalid files will be gracefully ignored by MultiQC."
        prefix: "A string for the MultiQC output directory: <prefix>/ and <prefix>.tar.gz"
        memory_gb: "RAM to allocate for task, specified in GB"
        disk_size_gb: "Disk space to allocate for task, specified in GB"
    }

    input {
        Array[File] input_files
        String prefix
        Int memory_gb = 5
        Int disk_size_gb = 20  # TODO why is this so high? Convert to dynamic
        Int max_retries = 1
    }

    String out_tar_gz = prefix + ".tar.gz"

    command <<<
        set -euo pipefail

        # set ENV variables for `multiqc`
        export LC_ALL=C.UTF-8
        export LANG=C.UTF-8

        echo "~{sep('\n', input_files)}" > file_list.txt

        multiqc -v \
            --no-ansi \
            -c multiqc_config.yaml \
            --file-list file_list.txt \
            -o ~{prefix}
        
        if [ ! -d ~{prefix} ]; then
            >&2 echo "MultiQC didn't find any valid files!"
            exit 1
        fi

        tar -czf ~{out_tar_gz} ~{prefix}
    >>>

    output {
        File multiqc_report = out_tar_gz
    }

    runtime {
        memory: "~{memory_gb} GB"
        disk: "~{disk_size_gb} GB"
        docker: 'quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0'
        maxRetries: max_retries
    }
}
