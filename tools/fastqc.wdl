## [Homepage](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

version 1.1

task fastqc {
    meta {
        description: "Generates a FastQC quality control metrics report for the input BAM file"
        outputs: {
            raw_data: "A zip archive of raw FastQC data. Can be parsed by MultiQC.",
            results: "A gzipped tar archive of all FastQC output files",
        }
    }

    parameter_meta {
        bam: "Input BAM format file to run FastQC on"
        prefix: "Prefix for the FastQC results directory. The extension `.tar.gz` will be added."
        use_all_cores: {
            description: "Use all cores? Recommended for cloud environments.",
            group: "common",
        }
        ncpu: {
            description: "Number of cores to allocate for task",
            group: "common",
        }
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        String prefix = basename(bam, ".bam") + ".fastqc_results"
        Boolean use_all_cores = false
        Int ncpu = 4
        Int modify_disk_size_gb = 0
    }

    String sanitized_prefix = sub(prefix, "[^a-zA-Z0-9_.-]", "_")
    String out_tar_gz = sanitized_prefix + ".tar.gz"
    String bam_basename = basename(bam, ".bam")
    String expected_zip_name = bam_basename + "_fastqc.zip"

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size * 2) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi

        mkdir ~{sanitized_prefix}
    
        fastqc -f bam \
            -o ~{sanitized_prefix} \
            -t "$n_cores" \
            ~{bam}
            
        # Create a checkpoint for the output exist
        if [ ! -f "~{sanitized_prefix}/~{expected_zip_name}" ]; then
            echo "ERROR: FastQC output file not found at expected location: ~{sanitized_prefix}/~{expected_zip_name}"
            # List contents to aid debugging
            ls -la ~{sanitized_prefix}/
            exit 1
        fi

        # tar achieve
        tar -czf ~{out_tar_gz} ~{sanitized_prefix}
    >>>

    output {
        File raw_data = "~{sanitized_prefix}/~{bam_basename}_fastqc.zip"
        File results = out_tar_gz
    }

    runtime {
        cpu: ncpu
        memory: "4 GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
        maxRetries: 1
    }
}


# Problem: if the prefic contain the special character or delimeter, it can casue the file path to interpreted incorectly
#Solution : we trying to use the prefix with sanitized_prefix --> update output section and validation

