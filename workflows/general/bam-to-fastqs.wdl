# SPDX-License-Identifier: MIT
# Copyright St. Jude Children's Research Hospital
version 1.1

import "../../tools/fq.wdl"
import "../../tools/samtools.wdl"

workflow bam_to_fastqs {
    meta {
        description: "Converts an input BAM file to one or more FASTQ files, performing QC checks along the way"
        outputs {
            # TODO rename read1s and read2s?
            read1s: "Array of FASTQ files corresponding to either `first` reads (if `paired_end = true`) or all reads (if `paired_end = false`)"
            read2s: "Array of FASTQ files corresponding to `last` reads (if `paired_end = true`)"
        }
    }

    parameter_meta {
        bam: "BAM file to split into FASTQs"
        paired_end: "Is the data paired-end (true) or single-end (false)?"
        use_all_cores: "Use all cores for multi-core steps?"
        max_retries: "Number of times to retry failed steps. Overrides task level defaults."
    }

    input {
        File bam
        Boolean paired_end = true
        Boolean use_all_cores = false
        Int? max_retries
    }

    call samtools.quickcheck { input: bam=bam, max_retries=max_retries }
    call samtools.split { input: bam=bam, use_all_cores=use_all_cores, max_retries=max_retries }
    scatter (split_bam in split.split_bams) {
        call samtools.collate_to_fastq as bam_to_fastq { input:
            bam=split_bam,
            paired_end=paired_end,
            interleaved=false,  # matches default but prevents user from overriding
            use_all_cores=use_all_cores,
            max_retries=max_retries
        }
    }

    scatter (reads in 
        zip(bam_to_fastq.read_one_fastq_gz, bam_to_fastq.read_two_fastq_gz)
    ) {
        call fq.fqlint { input:
            read_one_fastq=select_first([reads.left, "undefined"]),
            read_two_fastq=reads.right,
            max_retries=max_retries
        }
    }

    output {
        Array[File] read1s = select_all(bam_to_fastq.read_one_fastq_gz)
        Array[File?] read2s = bam_to_fastq.read_two_fastq_gz
    }
}
