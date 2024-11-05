version 1.1

import "../../data_structures/read_group.wdl"
import "../../tools/picard.wdl"
import "../general/bam-to-fastqs.wdl" as bam_to_fastqs_wf
import "./fastq-to-ubam-core.wdl" as ubam_core

workflow bam_to_ubam {
    meta {
        description: "Generate unaligned BAM file from aligned BAM file."
        outputs: {
            unaligned_bam: "Unaligned BAM file",
        }
        allowNestedInputs: true
    }

    parameter_meta {
        bam: "Input BAM from which to remove alignment information."
        prefix: "Prefix for the BAM file. The extension `.bam` will be added."
        validate_input: "Ensure input BAM is well-formed before beginning harmonization?"
        use_all_cores: "Use all cores? Recommended for cloud environments."
    }

    input {
        File bam
        File bwa_db
        File? restriction_sites
        String genome_id = "hg38"
        String prefix = basename(bam, ".bam")
        Boolean validate_input = true
        Boolean use_all_cores = false
    }

    if (validate_input) {
        call picard.validate_bam as validate_input_bam { input:
            bam,
        }
    }

    call read_group.get_read_groups { input:
        bam,
    }

    call bam_to_fastqs_wf.bam_to_fastqs { input:
        bam,
        paired_end = true,  # matches default but prevents user from overriding
        use_all_cores,
    }

    call ubam_core.fastq_to_ubam_core { input:
        read_one_fastqs_gz = bam_to_fastqs.read1s,
        read_two_fastqs_gz = select_all(bam_to_fastqs.read2s),
        read_groups = get_read_groups.read_groups,
        prefix,
        use_all_cores,
    }

    output {
        File unaligned_bam = fastq_to_ubam_core.unaligned_bam
    }
}
