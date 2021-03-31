version 1.0

import "https://raw.githubusercontent.com/stjudecloud/workflows/rnaseq-standard/v2.1.0/tools/estimate.wdl"

workflow ESTIMATE {
    input {
        File counts_file
        File gene_lengths_file
    }

    call estimate.calc_tpm { input: counts=counts_file, gene_lengths=gene_lengths_file }
    call estimate.run_ESTIMATE { input: gene_expression_file=calc_tpm.out }

    output {
        File tpm=calc_tpm.out
        File estimate_out=run_ESTIMATE.out
    }
}
