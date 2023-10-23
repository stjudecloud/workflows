## **WARNING:** this workflow is experimental! Use at your own risk!
#
# SPDX-License-Identifier: MIT
# Copyright St. Jude Children's Research Hospital
version 1.1

import "../../tools/bwa.wdl"
import "../../tools/picard.wdl"
import "../../tools/samtools.wdl"
import "../../tools/util.wdl"
import "../general/bam-to-fastqs.wdl" as bam_to_fastqs_wf

workflow dnaseq_standard_experimental {
    input {
        File bam
        File bwa_db
        Int? max_retries
        String prefix = basename(bam, ".bam")
        Boolean validate_input = true
        Boolean use_all_cores = false
    }

    if (validate_input) {
        call picard.validate_bam as validate_input_bam { input:
            bam=bam,
            max_retries=max_retries
        }
    }

    call util.get_read_groups { input:
        bam=bam,
        format_for_star=false,
        max_retries=max_retries
    }  # TODO what happens if no RG records?
    call bam_to_fastqs_wf.bam_to_fastqs { input:
        bam=bam,
        paired_end=true,  # matches default but prevents user from overriding
        use_all_cores=use_all_cores,
        max_retries=max_retries
    }

    scatter (tuple in zip(
        zip(bam_to_fastqs.read1s, bam_to_fastqs.read2s),
        get_read_groups.read_groups
    )) {
        call bwa.bwa_mem { input:
            read_one_fastq_gz=tuple.left.left,
            read_two_fastq_gz=tuple.left.right,
            bwa_db_tar_gz=bwa_db,
            # find tab literals, replace with '\\t' (which must be written as '\\\\t')
            # '\\t' is subbed into command blocks as '\t'
            read_group=sub(tuple.right, "\t", "\\\\t"),
            use_all_cores=use_all_cores,
            max_retries=max_retries
        }
        call picard.sort { input:
            bam=bwa_mem.bam,
            max_retries=max_retries
        }
    }
    call samtools.merge { input:
        bams=sort.sorted_bam,
        prefix=prefix,
        combine_pg=false,
        use_all_cores=use_all_cores,
        max_retries=max_retries
    }
}