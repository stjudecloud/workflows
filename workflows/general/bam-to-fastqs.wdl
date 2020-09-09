## # Bam to FastQs
##
## This WDL workflow converts an input BAM file to a set of fastq files for read 1 and read 2.
## It performs QC checks along the way to validate the input and output.
##
## ### Output:
##
## read1s
## : an array of files with the first read in the pair
##
## read2s
## : an array of files with the second read in the pair
##
## ## LICENSING:
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


import "https://raw.githubusercontent.com/stjudecloud/workflows/gzip-fastqs/tools/samtools.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/gzip-fastqs/tools/picard.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/gzip-fastqs/tools/fq.wdl"

workflow bam_to_fastqs {
    input {
        File bam
        Int samtools_sort_ncpu = 1
        Int bam_to_fastq_memory_gb = 40
        Int max_retries = 1
        Boolean detect_nproc = false
    }

    parameter_meta {
        bam: "BAM file to split into fastqs."
        samtools_sort_ncpu: "Number of CPUs to use while sorting the BAM."
        bam_to_fastq_memory_gb: "How much memory to provide while converting to fastqs."
        max_retries: "Maximum number of times to retry on a failure."
    }

    call samtools.quickcheck { input: bam=bam, max_retries=max_retries }
    call samtools.split { input: ncpu=samtools_sort_ncpu, bam=bam, max_retries=max_retries, detect_nproc=detect_nproc }
    scatter (split_bam in split.split_bams) {
        call picard.bam_to_fastq { input: bam=split_bam, memory_gb=bam_to_fastq_memory_gb, max_retries=max_retries }
    }
    scatter (reads in zip(bam_to_fastq.read1, bam_to_fastq.read2)) {
        call fq.fqlint { input: read1=reads.left, read2=reads.right, max_retries=max_retries }
    }

    output {
        Array[File] read1s = bam_to_fastq.read1
        Array[File] read2s = bam_to_fastq.read2
    }
}
