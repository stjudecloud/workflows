# SPDX-License-Identifier: MIT

version 1.1

import "../../tools/estimate.wdl"

workflow ESTIMATE {
    meta {
        description: "Runs the ESTIMATE software package on a feature counts file"
        external_help: "https://bioinformatics.mdanderson.org/estimate/"
        outputs: {
            tpm: "Transcripts Per Million file"
            estimate_result: "Final output of ESTIMATE"
        }
    }

    parameter_meta {
        counts_file: "A two column headerless TSV file with gene names in the first column and counts (as integers) in the second column. Entries starting with '__' will be discarded. Can be generated with `htseq.wdl`."
        gene_lengths_file: "A two column headered TSV file with gene names (matching those in the `counts` file) in the first column and feature lengths (as integers) in the second column. Can be generated with `calc-gene-lengths.wdl`."
        max_retries: "Number of times to retry failed steps. Overrides task level defaults."
    }

    input {
        File counts_file
        File gene_lengths_file
        Int? max_retries
    }

    call estimate.calc_tpm { input:
        counts=counts_file,
        gene_lengths=gene_lengths_file,
        max_retries=max_retries
    }
    call estimate.run_ESTIMATE { input:
        gene_expression_file=calc_tpm.tpm_file,
        max_retries=max_retries
    }

    output {
        File tpm=calc_tpm.tpm_file
        File estimate_result=run_ESTIMATE.estimate_file
    }
}
