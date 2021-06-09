## # RNASeq STAR DB build
##
## This WDL workflow generates a set of genome reference files usable by the STAR aligner from an input reference file in FASTA format.  
##
## ### Output
##
## reference_fa
## : the reference FASTA file
##
## gtf
## : the reference GTF file
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

import "https://raw.githubusercontent.com/stjudecloud/workflows/remove_read_string/tools/star.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/remove_read_string/tools/util.wdl"

workflow rnaseq_star_db_build {
    input {
        String reference_fa_url
        String gtf_url
        String reference_fa_name
        String gtf_name
    }

    parameter_meta {
        reference_fa_url: "URL to retrieve the reference FASTA file from"
        gtf_url: "URL to retrieve the reference GTF file from"
        reference_fa_name: "Name of output reference FASTA file"
        gtf_name: "Name of output GTF file"
    }

    call util.download as reference_download { input: url=reference_fa_url, outfilename=reference_fa_name }
    call util.download as gtf_download { input: url=gtf_url, outfilename=gtf_name }
    call star.build_db as star_db_build {
        input:
            reference_fasta=reference_download.outfile,
            gtf=gtf_download.outfile,
            stardb_dir_name="STARDB",
            ncpu=4,
    }

    output {
      File reference_fa = reference_download.outfile
      File gtf = gtf_download.outfile
      File stardb_tar_gz = star_db_build.stardb_out
    }
}
