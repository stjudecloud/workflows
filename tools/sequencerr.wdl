

version 1.0

task sequencerr {
    input {
        File bam
        File bai
        String prefix = basename(bam, ".bam")
        Boolean output_count_file = false
        Int max_retries = 1
    }

    parameter_meta {
        output_count_file: "Whether to output granular count file. Result file is large; uncompressed can reach over half input BAM size, compressed can be ~10% BAM size."
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil(bam_size * 2.2)

    String outfile = if output_count_file then prefix + ".sequencErr_results.tar.gz" else prefix + ".sequencErr.err"

    command {
        set -euo pipefail
        
        mv ~{bai} ~{bam}.bai || true
        sequencerr -pe=~{prefix}.sequencErr.err ~{bam} ~{prefix}.sequencErr.count

        if [ ~{output_count_file} == "true" ]; then
            mkdir ~{prefix}.sequencErr_results
            mv ~{prefix}.sequencErr.err ~{prefix}.sequencErr.count ~{prefix}.sequencErr_results
            tar -czf ~{prefix}.sequencErr_results.tar.gz ~{prefix}.sequencErr_results/
        fi
    }

    runtime {
        disk: disk_size + " GB"
        memory: "8 GB"
        docker: 'stjudecloud/sequencerr:branch-sequencErr-1.0.0'
        maxRetries: max_retries
    }

    output {
        File results = outfile
    }
}