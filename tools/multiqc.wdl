## [Homepage](https://multiqc.info/)

version 1.1

task multiqc {
    meta {
        description: "Generates a MultiQC quality control metrics report summary from input QC result files"
        outputs: {
            multiqc_report: "A gzipped tar archive of all MultiQC output files"
        }
    }

    parameter_meta {
        input_files: "An array of files for MultiQC to compile into a report. Invalid files will be gracefully ignored by MultiQC."
        prefix: "A string for the MultiQC output directory: <prefix>/ and <prefix>.tar.gz"
        config: "YAML file for configuring generated report"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        Array[File] input_files
        String prefix
        File? config
        Int modify_disk_size_gb = 0
    }

    Float input_size = size(input_files, "GiB")
    Int disk_size_gb = ceil(input_size) + 10 + modify_disk_size_gb
    String out_tar_gz = prefix + ".tar.gz"

    command <<<
        set -euo pipefail

        # set ENV variables for `multiqc`
        export LC_ALL=C.UTF-8
        export LANG=C.UTF-8

        echo "~{sep("\n", input_files)}" > file_list.txt

        # --strict is too strict. It causes errors due
        # to how our config adds 'custom-content' to the report.
        # Leaving this here as a warning not to try putting it back.
        # --require-logs might be useful at some point, but as of now,
        # it would cause errors. It could replace the check currently
        # run after multiqc is finished.
        # TODO: lots of other options to consider supporting.
        multiqc -v \
            --no-ansi \
            ~{if defined(config) then "-c " + config else ""} \
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
        memory: "4 GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/multiqc:1.22.2--pyhdfd78af_0"
        maxRetries: 1
    }
}
