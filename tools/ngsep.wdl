version 1.2

task germline_variant {
    meta {
        description: "Call germline variants using NGSEP"
        outputs: {
            vcf_output: "VCF file containing called germline variants"
        }
    }

    parameter_meta {
        reference_fasta: "Reference genome in FASTA format"
        bam: "Input BAM file with aligned reads"
        output_prefix: "Prefix for the output file with called variants"
        threads: "Number of threads to use"
        modify_disk_size_gb: "Additional disk size in GB to allocate"
    }

    input {
        File reference_fasta
        File bam
        String output_prefix = "ngsep_germline_output"
        Int threads = 4
        Int modify_disk_size_gb = 0
    }

    Int disk_size_gb = ceil(size(reference_fasta, "GiB") * 2)
        + ceil(size(bam, "GiB"))
        + 20
        + modify_disk_size_gb

    command <<<
        set -euo pipefail

        ref_fasta=~{basename(reference_fasta, ".gz")}
        gunzip -c "~{reference_fasta}" > "$ref_fasta" \
            || ln -sf "~{reference_fasta}" "$ref_fasta"

        java -Xmx16g -jar /usr/local/bin/NGSEPcore.jar \
            SingleSampleVariantsDetector \
            -r "$ref_fasta" \
            -i "~{bam}" \
            -o "~{output_prefix}" \
            -t "~{threads}"

        rm -rf "$ref_fasta"
    >>>

    output {
        Array[File] vcf_output = glob("~{output_prefix}*")
    }

    requirements {
        container: "ghcr.io/stjudecloud/ngsep:branch-minimap2-5.1.0-0"
        cpu: threads
        memory: "20 GB"
        disks: "~{disk_size_gb} GB"
    }
}
