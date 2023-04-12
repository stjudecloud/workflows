## # STAR DB build
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

import "../../tools/star.wdl"
import "../../tools/util.wdl"

workflow star_db_build {
    input {
        String reference_fa_url
        String reference_fa_name
        String? reference_fa_md5
        String gtf_url
        String gtf_name
        String? gtf_md5
        Int? max_retries
    }

    parameter_meta {
        reference_fa_url: "URL to retrieve the reference FASTA file from"
        reference_fa_name: "Name of output reference FASTA file"
        reference_fa_md5: "Expected md5sum of reference FASTA file"
        gtf_url: "URL to retrieve the reference GTF file from"
        gtf_name: "Name of output GTF file"
        gtf_md5: "Expected md5sum of GTF file"
        max_retries: "Number of times to retry failed steps. Overrides task level defaults."
    }

    call util.download as reference_download { input:
        url=reference_fa_url,
        outfile_name=reference_fa_name,
        md5sum=reference_fa_md5,
        max_retries=max_retries
    }
    call util.download as gtf_download { input:
        url=gtf_url,
        outfile_name=gtf_name,
        md5sum=gtf_md5,
        max_retries=max_retries
    }
    call star.build_star_db { input:
        reference_fasta=reference_download.outfile,
        gtf=gtf_download.outfile,
        max_retries=max_retries
    }

    output {
      File reference_fa = reference_download.outfile
      File gtf = gtf_download.outfile
      File stardb_tar_gz = build_star_db.stardb_out
    }
}
