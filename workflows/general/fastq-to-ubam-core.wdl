version 1.1

import "../../data_structures/read_group.wdl"
import "../../tools/picard.wdl"
import "../general/samtools-merge.wdl" as samtools_merge_wf

workflow fastq_to_ubam_core {
    meta {
        description: "Core pipeline to produce unaligned BAM from FASTQ files."
        outputs: {
            unaligned_bam: "Unaligned BAM file",
        }
        allowNestedInputs: true
    }

    parameter_meta {
        read_one_fastqs_gz: "Array of gzipped FASTQ files with 1st reads in pair. Assumes 1 file pair per read group."
        read_two_fastqs_gz: "Array of gzipped FASTQ files with 2nd reads in pair. Assumes 1 file pair per read group."
        read_groups: "An array of ReadGroup objects containing the read group information to output in the BAM file. One entry per file in `read_one_fastqs_gz`/`read_two_fastqs_gz`."
        prefix: "Prefix for the BAM file. The extension `.bam` will be added."
        validate_input: "Ensure input BAM is well-formed before beginning harmonization?"
        use_all_cores: "Use all cores? Recommended for cloud environments."
    }

    input {
        Array[File] read_one_fastqs_gz
        Array[File] read_two_fastqs_gz
        Array[ReadGroup] read_groups
        String prefix
        Boolean validate_input = true
        Boolean use_all_cores = false
    }

    scatter (tuple in zip(
        zip(read_one_fastqs_gz, read_two_fastqs_gz),
        read_groups
    )) {
        call picard.fastq_to_sam { input:
            read_one_fastq_gz = tuple.left.left,
            read_two_fastq_gz = tuple.left.right,
            prefix = basename(tuple.left.left, ".fastq.gz"),
            read_group_name = select_first([tuple.right.ID, ""]),
            sample_name = select_first([tuple.right.SM, ""]),
            library_name = select_first([tuple.right.LB, ""]),
            sequencing_center = select_first([tuple.right.CN, ""]),
            run_date = select_first([tuple.right.DT, ""]),
            platform_unit = select_first([tuple.right.PU, ""]),
            platform = select_first([tuple.right.PL, ""]),
            platform_model = select_first([tuple.right.PM, ""]),
        }
    }

    call samtools_merge_wf.samtools_merge as merge_unaligned_bam { input:
        bams = fastq_to_sam.unaligned_bam,
        prefix,
        use_all_cores,
    }

    output {
        File unaligned_bam = merge_unaligned_bam.merged_bam
    }
}
