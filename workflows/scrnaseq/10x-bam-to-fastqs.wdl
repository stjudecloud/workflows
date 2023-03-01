## # Cell Ranger Bam to FastQs
##
## This WDL workflow converts an input BAM file to a set of fastq files.
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
## fastqs
## : an array of files sufficient for localizing in Cell Ranger's expected format
##
## fastqs_archive
## : a compressed archive comtaining the array of fastq files
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


import "https://raw.githubusercontent.com/stjudecloud/workflows/master/tools/samtools.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/master/tools/cellranger.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/master/tools/fq.wdl"

workflow cell_ranger_bam_to_fastqs {
    input {
        File bam
        String pairing = "Paired-end"
        Boolean cellranger11 = false
        Boolean longranger20 = false
        Boolean gemcode = false
        Int bam_to_fastq_memory_gb = 40
        Boolean detect_nproc = false
        Int max_retries = 1
    }

    parameter_meta {
        bam: "BAM file to split into fastqs."
        bam_to_fastq_memory_gb: "How much memory to provide while converting to fastqs."
        max_retries: "Maximum number of times to retry on a failure."
    }

    call samtools.quickcheck { input: bam=bam, max_retries=max_retries }
    call cellranger.bamtofastq { input: bam=bam, cellranger11=cellranger11, longranger20=longranger20, gemcode=gemcode memory_gb=bam_to_fastq_memory_gb, max_retries=max_retries }
    scatter (reads in zip(bamtofastq.read1, bamtofastq.read2)) {
        if (pairing == "Paired-end") {
            call fq.fqlint as fqlint_pair { input: read1=reads.left, read2=reads.right, max_retries=max_retries }
        }
        if (pairing == "Single-end") {
            call fq.fqlint as fqlint_single { input: read1=reads.left, max_retries=max_retries }
        }
    }

    output {
        Array[File] fastqs = bamtofastq.fastqs
        File fastqs_archive = bamtofastq.fastqs_archive
        Array[File] read1s = bamtofastq.read1
        Array[File] read2s = bamtofastq.read2
    }
}

task parse_input {
    input {
        Boolean cellranger11
        Boolean longranger20
        Boolean gemcode
    }

    Int exclusive_arg = (if cellranger11 then 1 else 0) + (if longranger20 then 1 else 0) + (if gemcode then 1 else 0)

    command <<<
        if [ "~{exclusive_arg}" -gt 1 ]; then
            >&2 echo "Only one of cellranger11, longranger20, or gemcode can be set"
            exit 1
        fi
    >>>

    runtime {
        memory: "4 GB"
        disk: "1 GB"
        docker: 'ghcr.io/stjudecloud/util:1.2.0'
    }

    output {
        String input_check = "passed"
    }
}