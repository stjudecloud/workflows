## Description: 
##
## This WDL workflow generates a set of genome reference files usable by the STAR aligner from an input reference file in FASTA format.  
##
## Inputs: 
##
## reference_fasta - the genome for which to generate STAR reference files in FASTA format 
## gencode_gtf - the gene model file for the reference genome to use when generating STAR reference files 
## stardb_dir_name - a name to use for the output STAR reference directory
##
## Output: 
##
## out_dir - a directory containing the STAR reference files corresponding to the input genome 
##
## LICENSING :
## MIT License
##
## Copyright 2019 St. Jude Children's Research Hospital
##
## Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
##
##The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
##
##THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import "https://raw.githubusercontent.com/stjudecloud/workflows/jobin/rnaseq_v2_azure/tools/star.wdl"

workflow build_db {
    File reference_fasta
    File gencode_gtf 
    String stardb_dir_name

    call star.build_db {
        input:
            reference_fasta=reference_fasta,
            gencode_gtf=gencode_gtf,
            stardb_dir_name=stardb_dir_name,
            ncpu=4
    }

    output {
        # File out_dir = build_db.dir
        File stardb_zip = build_db.zip
    }
}
