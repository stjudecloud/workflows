## **WARNING:** this workflow is experimental! Use at your own risk!

version 1.1

import "../../data_structures/read_group.wdl"
import "../../tools/picard.wdl"
import "../../tools/samtools.wdl"
import "../general/bam-to-fastqs.wdl" as bam_to_fastqs_wf
import "./dnaseq-core.wdl" as dnaseq_core_wf

workflow dnaseq_standard_experimental {
    meta {
        name: "DNA-Seq Standard (Experimental)"
        description: "Aligns DNA reads using bwa"
        category: "Harmonization"
        outputs: {
            harmonized_bam: "Harmonized DNA-Seq BAM, aligned with bwa",
            harmonized_bam_index: "Index for the harmonized DNA-Seq BAM file",
        }
        allowNestedInputs: true
    }

    parameter_meta {
        bam: "Input BAM to realign"
        bwa_db: "Gzipped tar archive of the bwa reference files. Files should be at the root of the archive."
        sample_override: "Value to override the SM field of *every* read group."
        prefix: "Prefix for the BAM file. The extension `.bam` will be added."
        aligner: {
            description: "BWA aligner to use",
            choices: [
                "mem",
                "aln"
            ],
        }
        validate_input: "Ensure input BAM is well-formed before beginning harmonization?"
        use_all_cores: "Use all cores? Recommended for cloud environments."
        reads_per_file: "Controls the number of reads per FASTQ file for internal split to run BWA in parallel."
        subsample_n_reads: "Only process a random sampling of `n` reads. Any `n`<=`0` for processing entire input."
    }

    input {
        File bam
        File bwa_db
        String? sample_override
        String prefix = basename(bam, ".bam")
        String aligner = "mem"
        Boolean validate_input = true
        Boolean use_all_cores = false
        Int reads_per_file = 10000000
        Int subsample_n_reads = -1
    }

    call parse_input { input:
        aligner
    }

    if (validate_input) {
        call picard.validate_bam as validate_input_bam { input:
            bam,
        }
    }

    if (subsample_n_reads > 0) {
        call samtools.subsample after parse_input { input:
            bam,
            desired_reads = subsample_n_reads,
            use_all_cores,
        }
    }
    File selected_bam = select_first([subsample.sampled_bam, bam])

    call read_group.get_read_groups { input:
        bam = selected_bam,
    }

    call bam_to_fastqs_wf.bam_to_fastqs { input:
        bam = selected_bam,
        paired_end = true,  # matches default but prevents user from overriding
        use_all_cores,
    }

    call dnaseq_core_wf.dnaseq_core_experimental { input:
        read_one_fastqs_gz = bam_to_fastqs.read1s,
        read_two_fastqs_gz = select_all(bam_to_fastqs.read2s),
        bwa_db,
        reads_per_file,
        read_groups = get_read_groups.read_groups,
        prefix,
        aligner,
        use_all_cores,
        sample_override,
    }

    output {
        File harmonized_bam = dnaseq_core_experimental.harmonized_bam
        File harmonized_bam_index = dnaseq_core_experimental.harmonized_bam_index
    }
}

task parse_input {
    meta {
        description: "Parses and validates the `dnaseq_standard` workflow's provided inputs"
        outputs: {
            check: "Dummy output to indicate success and to enable call-caching"
        }
    }

    parameter_meta {
        aligner: {
            description: "BWA aligner to use",
            choices: [
                "mem",
                "aln"
            ],
        }
    }

    input {
        String aligner
    }

    command <<<
        if [ "~{aligner}" != "mem" ] \
            && [ "~{aligner}" != "aln" ]
        then
            >&2 echo "Aligner must be:"
            >&2 echo "'mem' or 'aln'"
            exit 1
        fi
    >>>

    output {
        String check = "passed"
    }

    runtime {
        memory: "4 GB"
        disks: "10 GB"
        container: "ghcr.io/stjudecloud/util:branch-docker-1.5.0"
        maxRetries: 0
    }
}
