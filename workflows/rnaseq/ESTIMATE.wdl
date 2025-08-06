version 1.1

import "../../tools/estimate.wdl"
import "../../tools/htseq.wdl"

workflow estimate {
    meta {
        name: "ESTIMATE"
        description: "Runs the ESTIMATE software package on a feature counts file"
        warning: "**[DEPRECATED]**"
        external_help: "https://bioinformatics.mdanderson.org/estimate/"
        outputs: {
            tpm: "Transcripts Per Million file",
            estimate_result: "Final output of ESTIMATE",
        }
        deprecated: true
    }

    parameter_meta {
        counts_file: "A two column headerless TSV file with gene names in the first column and counts (as integers) in the second column. Entries starting with '__' will be discarded. Can be generated with `htseq.wdl`."
        gene_lengths_file: "A two column headered TSV file with gene names (matching those in the `counts` file) in the first column and feature lengths (as integers) in the second column. Can be generated with `calc-gene-lengths.wdl`."
    }

    input {
        File counts_file
        File gene_lengths_file
    }

    call htseq.calc_tpm { input:
        counts = counts_file,
        gene_lengths = gene_lengths_file,
    }
    call estimate.run_estimate { input:
        gene_expression_file = calc_tpm.tpm_file,
    }

    output {
        File tpm = calc_tpm.tpm_file
        File estimate_result = run_estimate.estimate_file
    }
}
