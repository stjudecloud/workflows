## # Make QC Reference
##
## Download/create all the reference files needed for `quality-check-standard.wdl`.
## This includes: a reference FASTA, a GTF, the database used by FastQ Screen, and exonic/coding regions BEDs
## for use in restricting mosdepth coverage analysis to selected regions.
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

import "../../tools/util.wdl"
import "../../tools/fastq_screen.wdl"
import "../../tools/mosdepth.wdl"

workflow make_qc_reference {
    input {
        String reference_fa_url
        String gtf_url
        String reference_fa_name
        String gtf_name
        Int max_retries = 1
    }

    parameter_meta {
        reference_fa_url: "URL to retrieve the reference FASTA file from"
        gtf_url: "URL to retrieve the reference GTF file from"
        reference_fa_name: "Name of output reference FASTA file"
        gtf_name: "Name of output GTF file"
    }

    call util.download as reference_download {
        input:
            url=reference_fa_url,
            outfilename=reference_fa_name,
            max_retries=max_retries
    }
    call util.download as gtf_download {
        input:
            url=gtf_url,
            outfilename=gtf_name,
            max_retries=max_retries
    }

    call fastq_screen.build_db as fastq_screen_build_db
    call mosdepth.make_coverage_regions_beds {
        input:
            gtf=gtf_download.outfile,
            max_retries=max_retries
    }

    output {
        File reference_fa = reference_download.outfile
        File gtf = gtf_download.outfile
        File fastq_screen_db = fastq_screen_build_db.db
        File exon_bed = make_coverage_regions_beds.exon_bed
        File CDS_bed = make_coverage_regions_beds.CDS_bed
    }
}