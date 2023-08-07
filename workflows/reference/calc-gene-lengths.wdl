## # Calculate Gene Lengths
##
## Calculates gene lengths from a GTF file.
## TODO should explain algorithm somewhere. Lots of different ways to calculate gene lengths
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

version 1.1

import "../../tools/util.wdl"

workflow calculate_gene_lengths {
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
        File gene_lengths=calc_gene_lens.gene_lengths
    }
}
