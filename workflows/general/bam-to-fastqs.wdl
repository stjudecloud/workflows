version 1.1

import "../../tools/fq.wdl"
import "../../tools/samtools.wdl"

workflow bam_to_fastqs {
    meta {
        description: "Converts an input BAM file to one or more FASTQ files, performing QC checks along the way"
        outputs: {
            read1s: "Array of FASTQ files corresponding to either `first` reads (if `paired_end = true`) or all reads (if `paired_end = false`)",
            read2s: "Array of FASTQ files corresponding to `last` reads (if `paired_end = true`)"
        }
        allowNestedInputs: true
    }

    parameter_meta {
        bam: "BAM file to split into FASTQs"
        paired_end: "Is the data Paired-End (true) or Single-End (false)?"
        use_all_cores: "Use all cores for multi-core steps?"
    }

    input {
        File bam
        Boolean paired_end = true
        Boolean use_all_cores = false
    }

    call samtools.quickcheck { input: bam = bam }
    call samtools.split { input:
        bam,
        use_all_cores,
    }
    scatter (split_bam in split.split_bams) {
        call samtools.bam_to_fastq { input:
            bam = split_bam,
            paired_end,
            interleaved = false,  # matches default but prevents user from overriding
            use_all_cores,
        }
    }
    scatter (reads in zip(bam_to_fastq.read_one_fastq_gz, bam_to_fastq.read_two_fastq_gz)) {
        call fq.fqlint { input:
            read_one_fastq = select_first([reads.left, "undefined"]),
            read_two_fastq = reads.right,
        }
    }

    output {
        Array[File] read1s = select_all(bam_to_fastq.read_one_fastq_gz)
        Array[File?] read2s = bam_to_fastq.read_two_fastq_gz
    }
}
