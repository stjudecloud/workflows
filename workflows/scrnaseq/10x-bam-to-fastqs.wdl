## # Cell Ranger BAM to FASTQs
##
## This WDL workflow converts an input BAM file to a set of FASTQ files.
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
## : a compressed archive containing the array of FASTQ files
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

version 1.1

import "../../tools/cellranger.wdl"
import "../../tools/fq.wdl"
import "../../tools/samtools.wdl"

workflow cell_ranger_bam_to_fastqs {
    meta {
        description: "Convert a 10x Genomics BAM file to FASTQs."
        allowNestedInputs: true
        outputs: {
            fastqs: "FASTQ files with reads.",
            fastqs_archive: "Compressed archive of FASTQ files.",
            read1s: "Gzipped read 1 FASTQ files.",
            read2s: "Gzipped read 2 FASTQ files.",
        }
    }
    parameter_meta {
        bam: "BAM file to split into FASTQs."
        cellranger11: "Convert a BAM produced by Cell Ranger 1.0-1.1"
        longranger20: "Convert a BAM produced by Longranger 2.0"
        gemcode: "Convert a BAM produced from GemCode data (Longranger 1.0 - 1.3)"
        use_all_cores: "Use all cores for multi-core steps?"
    }

    input {
        File bam
        Boolean cellranger11 = false
        Boolean longranger20 = false
        Boolean gemcode = false
        Boolean use_all_cores = false
    }

    call samtools.quickcheck { input: bam }
    call cellranger.bamtofastq { input:
        bam,
        cellranger11,
        longranger20,
        gemcode,
        use_all_cores,
    }
    scatter (reads in zip(bamtofastq.read_one_fastq_gz, bamtofastq.read_two_fastq_gz)) {
        call fq.fqlint { input:
            read_one_fastq = reads.left,
            read_two_fastq = reads.right,
        }
    }

    output {
        Array[File] fastqs = bamtofastq.fastqs
        File fastqs_archive = bamtofastq.fastqs_archive
        Array[File] read1s = bamtofastq.read_one_fastq_gz
        Array[File] read2s = bamtofastq.read_two_fastq_gz
    }
}

task parse_input {
    meta {
        description: "Parse 10x-bam-to-fastqs workflow inputs and validate"
    }
    input {
        Boolean cellranger11
        Boolean longranger20
        Boolean gemcode
    }
    parameter_meta {
        cellranger11: "Convert a BAM produced by Cell Ranger 1.0-1.1"
        longranger20: "Convert a BAM produced by Longranger 2.0"
        gemcode: "Convert a BAM produced from GemCode data (Longranger 1.0 - 1.3)"
    }

    Int exclusive_arg = (if cellranger11 then 1 else 0)
        + (if longranger20 then 1 else 0)
        + (if gemcode then 1 else 0)

    command <<<
        if [ "~{exclusive_arg}" -gt 1 ]; then
            >&2 echo "Only one of cellranger11, longranger20, or gemcode can be set"
            exit 1
        fi
    >>>

    output {
        String input_check = "passed"
    }

    runtime {
        memory: "4 GB"
        disk: "10 GB"
        container: "ghcr.io/stjudecloud/util:1.3.0"
        maxRetries: 1
    }
}
