## Description: 
##
## This WDL workflow runs the STAR aligner on a set of input fastq files.  
##
## Inputs: 
##
## read1s - an array of files with the first read in the pair
## read2s - an array of files with the second read in the pair
## stardb_dir - a STAR reference directory that was generated using STAR genomeGenerate
## output_prefix - the prefix name that should be used to name the STAR outputs
## read_groups (optional) - a list of the read group entries to use in the output BAM from STAR
##
## Output: 
##
## star_bam - the aligned reads in BAM format from STAR
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

workflow star_alignment {
    Array[File] read_one_fastqs
    Array[File] read_two_fastqs
    File stardb_zip
    String output_prefix
    String? read_groups

    call star.alignment {
        input:
            read_one_fastqs=read_one_fastqs,
            read_two_fastqs=read_two_fastqs,
            stardb_zip=stardb_zip,
            output_prefix=output_prefix,
            read_groups=read_groups
    }

    output {
       File star_bam = alignment.star_bam
    }
}
