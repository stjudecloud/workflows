# SPDX-License-Identifier: MIT

# TODO should explain algorithm somewhere. Lots of different ways to calculate gene lengths

version 1.1

import "../../tools/util.wdl"

workflow calculate_gene_lengths {
    meta {
        description: "Calculate gene lengths from a GTF feature file"  # TODO explain algorithm
        outputs: {
            gene_lengths: "A two column headered TSV file with gene names in the first column and feature lengths (as integers) in the second column"
        }
    }

    parameter_meta {
        gtf: "GTF feature file"
        max_retries: "Number of times to retry failed steps. Overrides task level defaults."
    }

    input {
        File gtf
        Int? max_retries
    }

    call util.calc_gene_lengths { input:
        gtf=gtf,
        max_retries=max_retries
    }

    output {
        File gene_lengths=calc_gene_lengths.gene_lengths
    }
}
