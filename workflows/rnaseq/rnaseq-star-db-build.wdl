## # RNASeq STAR DB build
##
## This WDL workflow generates a set of genome reference files usable by the STAR aligner from an input reference file in FASTA format.  
##
## ### Output
##
## reference_fa
## : the reference FASTA file
##
## gencode_gtf
## : the reference gencode GTF file
##
## stardb_tar_gz
## : the STAR DB folder in .tar.gz format
##
## ## LICENSING
##
## #### MIT License
##
## Copyright 2019 St. Jude Children's Research Hospital
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

import "https://raw.githubusercontent.com/stjudecloud/workflows/master/tools/star.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/master/tools/gzip.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/master/tools/wget.wdl"

workflow rnaseq_star_db_build {
    input {
        String reference_fa_url
        String gencode_gtf_url
    }

    parameter_meta {
        reference_fa_url: "URL to retrieve the reference FASTA file from."
        gencode_gtf_url: "URL to retrieve the reference gencode GTF file."
    }

    call wget.download as reference_download { input: url=reference_fa_url, outfilename="GRCh38_no_alt.fa.gz" }
    call gzip.unzip as reference_unzip { input: infile=reference_download.outfile }
    call wget.download as gencode_download { input: url=gencode_gtf_url, outfilename="gencode.v31.gtf.gz" }
    call gzip.unzip as gencode_unzip { input: infile=gencode_download.outfile }
    call star.build_db as star_db_build {
        input:
            reference_fasta=reference_unzip.outfile,
            gencode_gtf=gencode_unzip.outfile,
            stardb_dir_name="STARDB",
            ncpu=4,
    }

    output {
      File reference_fa = reference_unzip.outfile
      File gencode_gtf = gencode_unzip.outfile
      File stardb_tar_gz = star_db_build.stardb_out
    }
}
