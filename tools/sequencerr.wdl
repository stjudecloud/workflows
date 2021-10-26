## # sequencErr
##
## This WDL tool wraps the sequencErr software (unpublished).
## sequencErr approximates error rates in Illumina sequencing machines down to tile level precision.

version 1.0

task sequencerr {
    input {
        File bam
        File bai
        String? prefix = basename(bam, ".bam") + ".sequencErr"
        Boolean output_count_file = false
        Int max_retries = 1
    }

    parameter_meta {
        output_count_file: "Whether to output granular count file. Result file is large; uncompressed can reach over half input BAM size, compressed can be ~10% BAM size."
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil(bam_size * 2.2)

    # There's an odd bug in Cromwell requiring this convoluted structure
    String parsed_prefix = select_first([prefix, basename(bam, ".bam") + ".sequencErr"])
    String outfile = if output_count_file then parsed_prefix + "_results.tar.gz" else parsed_prefix + ".err"

    command {
        set -euo pipefail
        
        mv ~{bai} ~{bam}.bai || true
        sequencerr -pe=~{parsed_prefix}.err ~{bam} ~{parsed_prefix}.count

        if [ ~{output_count_file} == "true" ]; then
            mkdir ~{parsed_prefix}_results
            mv ~{parsed_prefix}.err ~{parsed_prefix}.count ~{parsed_prefix}_results
            tar -czf ~{parsed_prefix}_results.tar.gz ~{parsed_prefix}_results/
            rm -r ~{parsed_prefix}_results/
        else
            rm ~{parsed_prefix}.count
        fi
    }

    runtime {
        disk: disk_size + " GB"
        memory: "8 GB"
        docker: 'ghcr.io/stjudecloud/sequencerr:1.0.0'
        maxRetries: max_retries
    }

    output {
        File results = outfile
    }
}