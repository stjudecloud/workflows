## # MultiQC
##
## This WDL tool wraps the [MultiQC](https://multiqc.info/) tool.
## MultiQC aggregates quality control results for bioinformatics.

version 1.0

task multiqc {
    input {
        Array[File] input_files
        String output_prefix
        Array[String] extra_fn_clean_exts = []
        Int max_retries = 1
        Int memory_gb = 5
        Int disk_size = 20
    }

    String out_directory = output_prefix + ".multiqc"
    String out_tar_gz = out_directory + ".tar.gz"

    command {
        set -eo pipefail
        
        # set ENV variables for `multiqc`
        export LC_ALL=C.UTF-8
        export LANG=C.UTF-8
        
        echo "~{sep="\n" input_files}" > file_list.txt

        echo "~{sep="\n" extra_fn_clean_exts}" > extensions.txt

        echo "extra_fn_clean_exts:" > multiqc_config.yaml
        while read -r ext; do
            echo "    - $ext" >> multiqc_config.yaml
        done < extensions.txt

        multiqc -vvv -c multiqc_config.yaml \
            --file-list file_list.txt -o ~{out_directory}
        
        if [ ! -d ~{out_directory} ]; then
            >&2 echo "MultiQC didn't find any valid files!"
            exit 1
        fi

        tar -czf ~{out_tar_gz} ~{out_directory}
    }

    runtime {
        disk: disk_size + " GB"
        docker: 'quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0'
        memory: memory_gb + " GB"
        maxRetries: max_retries
    }

    output {
        File out = out_tar_gz
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "This WDL tool generates a MultiQC quality control metrics report summary from input QC result files."
    }

    parameter_meta {
        input_files: "A non-empty array of files for MultiQC to compile into a report. Invalid files will be gracefully ignored by MultiQC."
        output_prefix: "A string for the MultiQC output directory: <prefix>.multiqc/ and <prefix>.multiqc.tar.gz"
    }
}
