version 1.1

task bam2pairs {
    meta {
        description: "Convert BAM file to pairs file with Pairix"
        outputs: {
            pairs: "Mapped read pairs from BAM file"
        }
    }

    parameter_meta {
        bam: "BAM file from which to extra read pairs"
        prefix: "Prefix for the pairs file. The extension `.pairs.txt` will be added."
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        String prefix = basename(bam, ".bam")
        Int modify_disk_size_gb = 0
    }

    Float input_bam_size = size(bam, "GiB")
    Int disk_size_gb = (
        ceil((input_bam_size) * 2) + 10 + modify_disk_size_gb
    )

    command <<<
        bam2pairs \
            ~{bam} \
            ~{prefix}.pairs.txt
    >>>

    output {
        pairs = ~{prefix}.pairs.txt
    }

    runtime {
        cpu: ncpu
        memory: "4 GB"
        disk: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/pairix:0.3.7--py38h7cf9df2_5"
        maxRetries: 1
    }
}