version 1.1

import "../../tools/md5sum.wdl"
import "../../tools/picard.wdl"
import "../../tools/samtools.wdl"
import "https://raw.githubusercontent.com/stjude/XenoCP/4.0.0-alpha/wdl/workflows/xenocp.wdl" as xenocp_wf

workflow alignment_post {
    meta {
        description: "Runs a series of standard processing tools that should immediately follow alignment, regardless of data-type"
        outputs: {
            processed_bam: "Input BAM after being transformed by standard processing",
            bam_index: "BAI index associated with `processed_bam`",
            bam_checksum:  "STDOUT of the `md5sum` command run on the input BAM that has been redirected to a file"
        }
        allowNestedInputs: true
    }

    parameter_meta {
        bam: "Input BAM format file to process"
        mark_duplicates: "Add SAM flag to computationally determined duplicate reads?"
        contaminant_db: "A compressed reference database corresponding to the aligner chosen with `xenocp_aligner` for the contaminant genome"
        xenocp_aligner: {
            description: "Aligner to use to map reads to the host genome for detecting contamination",
            choices: [
                'bwa aln',
                'bwa mem',
                'star'
            ]
        }
        cleanse_xenograft: "If true, use XenoCP to unmap reads from contaminant genome"
        use_all_cores: "Use all cores for multi-core steps?"
    }

    input {
        File bam
        Boolean mark_duplicates
        File? contaminant_db
        String xenocp_aligner = ""
        Boolean cleanse_xenograft = false
        Boolean use_all_cores = false
    }

    call picard.sort as picard_sort { input: bam=bam }

    if (cleanse_xenograft) {
        call samtools.index as pre_xenocp_index { input:
            bam=picard_sort.sorted_bam,
            use_all_cores=use_all_cores,
        }

        call xenocp_wf.xenocp { input:
            input_bam=picard_sort.sorted_bam,
            input_bai=pre_xenocp_index.bam_index,
            reference_tar_gz=select_first([contaminant_db, ""]),
            aligner=xenocp_aligner,
            skip_duplicate_marking=true
        }
    }
    if (mark_duplicates) {
        call picard.mark_duplicates as picard_markdup { input:
            bam=select_first([xenocp.bam, picard_sort.sorted_bam]),
        }
    }

    File aligned_bam = select_first([
        picard_markdup.duplicate_marked_bam,
        xenocp.bam,
        picard_sort.sorted_bam
    ])

    call samtools.index as samtools_index { input:
        bam=aligned_bam,
        use_all_cores=use_all_cores,
    }
    File aligned_bam_index = samtools_index.bam_index

    # TODO: should we output the validater report?
    call picard.validate_bam { input: bam=aligned_bam }

    call md5sum.compute_checksum { input: file=aligned_bam }

    output {
        File processed_bam = aligned_bam
        File bam_index = aligned_bam_index
        File bam_checksum = compute_checksum.md5sum
    }
}
