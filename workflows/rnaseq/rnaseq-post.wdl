version 1.0

import "../../tools/htseq.wdl"
import "../../tools/ngsderive.wdl"

workflow rnaseq_post {
    input {
        File bam
        File bam_index
        File gtf
        String strandedness = ""
        Int max_retries = 1
    }

    parameter_meta {
        bam: "Input BAM format file"
        strandedness: "empty, 'Stranded-Reverse', 'Stranded-Forward', or 'Unstranded'. If missing, will be inferred"
        max_retries: "Number of times to retry failed steps"
    }

    call ngsderive.infer_strandedness as ngsderive_strandedness {
        input:
            bam=bam,
            bai=bam_index,
            gtf=gtf,
            max_retries=max_retries
    }
    String parsed_strandedness = read_string(ngsderive_strandedness.strandedness)

    call htseq.count as htseq_count { input:
        bam=bam,
        gtf=gtf,
        provided_strandedness=strandedness,
        inferred_strandedness=parsed_strandedness,
        max_retries=max_retries
    }

    output {
        File gene_counts = htseq_count.out
        File inferred_strandedness = ngsderive_strandedness.strandedness_file
    }
}