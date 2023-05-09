## # Cell Ranger Bam to FastQs
##
## This WDL workflow converts an input BAM file to a set of FastQ files.
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
## : a compressed archive comtaining the array of FastQ files
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


import "../../tools/samtools.wdl"
import "../../tools/cellranger.wdl"
import "../../tools/fq.wdl"

workflow cell_ranger_bam_to_fastqs {
    input {
        File bam
        Boolean paired = true
        Boolean cellranger11 = false
        Boolean longranger20 = false
        Boolean gemcode = false
        Int bam_to_fastq_memory_gb = 40
        Boolean use_all_cores = false
        Int? max_retries
    }

    parameter_meta {
        bam: "BAM file to split into FastQs."
        paired: "Is the data paired-end (true) or single-end (false)?"
        cellranger11: "Convert a BAM produced by Cell Ranger 1.0-1.1"
        longranger20: "Convert a BAM produced by Longranger 2.0"
        gemcode: "Convert a BAM produced from GemCode data (Longranger 1.0 - 1.3)"
        bam_to_fastq_memory_gb: "How much memory to provide while converting to FastQs."
        use_all_cores: "Use all cores for multi-core steps?"
        max_retries: "Number of times to retry failed steps. Overrides task level defaults."
    }

    call samtools.quickcheck { input: bam=bam, max_retries=max_retries }
    call cellranger.bamtofastq { input: bam=bam, cellranger11=cellranger11, longranger20=longranger20, gemcode=gemcode, memory_gb=bam_to_fastq_memory_gb, max_retries=max_retries }
    scatter (reads in zip(bamtofastq.read1, bamtofastq.read2)) {
        if (paired) {
            call fq.fqlint as fqlint_pair { input: read1=reads.left, read2=reads.right, max_retries=max_retries }
        }
        if (! paired) {
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
