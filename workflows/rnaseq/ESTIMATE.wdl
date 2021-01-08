version 1.0

import "https://raw.githubusercontent.com/stjudecloud/workflows/ESTIMATE/tools/estimate.wdl"

workflow ESTIMATE {
    input {
        File counts_file
        File gene_lengths_file
    }

    call estimate.calc_tpm { input: counts=counts_file, gene_lengths=gene_lengths_file }
    call estimate.run_ESTIMATE { input: gene_expression_file=calc_tpm.out }

    output {
        File gene_lengths=calc_tpm.out
        File estimate_unfiltered_out=run_ESTIMATE.unfiltered_out
        File estimate_filtered_out=run_ESTIMATE.filtered_out
    }
}