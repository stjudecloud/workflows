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
        Array[String] mosdepth_labels = []
        Int disk_size = 20
        Int memory_gb = 5
        Int max_retries = 1
    }

    parameter_meta {
        input_files: "An array of files for MultiQC to compile into a report. Invalid files will be gracefully ignored by MultiQC."
        output_prefix: "A string for the MultiQC output directory: <prefix>.multiqc/ and <prefix>.multiqc.tar.gz"
    }

    String out_directory = output_prefix + ".multiqc"
    String out_tar_gz = out_directory + ".tar.gz"

    command {
        set -euo pipefail
        
        # set ENV variables for `multiqc`
        export LC_ALL=C.UTF-8
        export LANG=C.UTF-8
        
        echo "~{sep="\n" input_files}" > file_list.txt

        # Start YAML generation
        echo "extra_fn_clean_exts:" > multiqc_config.yaml
        
        # if extra extensions are to be cleaned, add them to YAML
        if [ "~{if (length(extra_fn_clean_exts) > 0) then "true" else ""}" = "true" ]; then
            echo "~{sep="\n" extra_fn_clean_exts}" > extensions.txt
            while read -r ext; do
                echo "  - $ext"
            done < extensions.txt >> multiqc_config.yaml
        fi

        # if mosdepth labels are provided, add them to `extra_fn_clean_exts`
        if [ "~{if (length(mosdepth_labels) > 0) then "true" else ""}" = "true" ]; then
            # add default mosdepth label ("whole_genome") to 'labels.txt'
            echo "whole_genome" > labels.txt

            # add the rest of the provided mosdepth labels to 'labels.txt'
            echo "~{sep="\n" mosdepth_labels}" >> labels.txt
            while read -r label; do
                echo "  - .$label"
            done < labels.txt >> multiqc_config.yaml

            # next YAML section
            echo "top_modules:" >> multiqc_config.yaml
            
            # create a "top module" entry for each mosdepth label provided
            while read -r label; do
                echo "  - mosdepth:"
                echo "      name: \"mosdepth ($label)\""
                echo "      path_filters:"
                echo "        - \"*.$label.mosdepth.*\""
            done < labels.txt >> multiqc_config.yaml
        fi

        echo "*** Generated Config ***"
        cat multiqc_config.yaml
        echo "*** End ***"

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
}
