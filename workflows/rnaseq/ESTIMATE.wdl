## # ESTIMATE
##
## Runs the [ESTIMATE software package](https://bioinformatics.mdanderson.org/estimate/) on a feature counts file.
##
## ## LICENSING:
## 
## #### MIT License
##
## Copyright 2020-Present St. Jude Children's Research Hospital
##
## Permission is hereby granted, free of charge, to any person obtaining a copy of this
## software and associated documentation files (the "Software"), to deal in the Software
## without restriction, including without limitation the rights to use, copy, modify, merge,
## publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons
## to whom the Software is furnished to do so, subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in all copies or
## substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
## BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
## NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
## DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
## OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

version 1.0

import "../../tools/estimate.wdl"

workflow ESTIMATE {
    input {
        File counts_file
        File gene_lengths_file
        Int? max_retries
    }

    parameter_meta {
        counts: "A two column headerless TSV file with gene names in the first column and counts (as integers) in the second column. Entries starting with '__' will be discarded. Can be generated with `htseq.wdl`."
        gene_lengths: "A two column headered TSV file with gene names (matching those in the `counts` file) in the first column and feature lengths (as integers) in the second column. Can be generated with `calc-gene-lengths.wdl`."
        max_retries: "Number of times to retry failed steps. Overrides task level defaults."
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
        File estimate_out=run_ESTIMATE.estimate_file
    }
}
