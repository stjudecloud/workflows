## # BWA DB build
##
## This WDL workflow generates a set of genome reference files usable by the BWA aligner from an input reference file in FASTA format.  
##
## ### Output
##
## reference_fa
## : the reference FASTA file
##
## bwadb_tar_gz
## : the BWA reference folder in .tar.gz format
##
## ## LICENSING
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

import "../../tools/bwa.wdl"
import "../../tools/util.wdl"

workflow bwa_db_build {
    parameter_meta {
        reference_fa_url: "URL to retrieve the reference FASTA file from."
        reference_fa_name: "Name of output reference FASTA file"
        reference_fa_md5: "Expected md5sum of reference FASTA file"
        max_retries: "Number of times to retry failed steps"
    }

    input {
        String reference_fa_url
        String reference_fa_name
        String? reference_fa_md5
        Int max_retries = 1
    }

    call util.download as reference_download {
        input:
            url=reference_fa_url,
            outfile_name=reference_fa_name,
            md5sum=reference_fa_md5,
            max_retries=max_retries
    }
    call bwa.build_bwa_db {
        input:
            reference_fasta=reference_download.outfile,
            max_retries=max_retries
    }

    output {
      File reference_fa = reference_download.outfile
      File bwadb_tar_gz = build_bwa_db.bwadb_tar_gz
    }
}
