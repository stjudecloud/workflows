## **WARNING:** this workflow is experimental! Use at your own risk!

version 1.1

import "../../data_structures/read_group.wdl"
import "../../tools/picard.wdl"
import "../../tools/util.wdl"
import "../general/bam-to-fastqs.wdl" as bam_to_fastqs_wf
import "./dnaseq-core.wdl" as dnaseq_core_wf

workflow dnaseq_standard_experimental {
    meta {
        description: "Aligns DNA reads using bwa"
        outputs: {
            harmonized_bam: "Harmonized DNA-Seq BAM, aligned with bwa"
            harmonized_bam_index: "Index for the harmonized DNA-Seq BAM file"
        }
        allowNestedInputs: true
    }
    parameter_meta {
        bam: "Input BAM to realign"
        bwa_db: "Gzipped tar archive of the bwa reference files. Files should be at the root of the archive."
        reads_per_file: "Number of reads per FASTQ file to output."
        prefix: "Prefix for the BAM file. The extension `.bam` will be added."
        aligner: {
            description: "BWA aligner to use",
            choices: ["mem", "aln"]
        }
        validate_input: "Ensure input BAM is well-formed before beginning harmonization?"
        use_all_cores: "Use all cores? Recommended for cloud environments."
    }
    input {
        File bam
        File bwa_db
        Int reads_per_file = 10000000
        String prefix = basename(bam, ".bam")
        String aligner = "mem"
        Boolean validate_input = true
        Boolean use_all_cores = false
    }

    call parse_input { input:
        aligner=aligner
    }

    if (validate_input) {
        call picard.validate_bam as validate_input_bam { input:
            bam=bam,
        }
    }

    # call util.get_read_groups { input:
    #     bam=bam,
    #     format_for_star=false,
    # }  # TODO what happens if no RG records?
    call read_group.get_ReadGroups { input:
        bam=bam,
    }

    call bam_to_fastqs_wf.bam_to_fastqs { input:
        bam=bam,
        paired_end=true,  # matches default but prevents user from overriding
        use_all_cores=use_all_cores,
    }

    call dnaseq_core_wf.dnaseq_core_experimental { input:
        read_one_fastqs_gz = bam_to_fastqs.read1s,
        read_two_fastqs_gz = select_all(bam_to_fastqs.read2s),
        bwa_db = bwa_db,
        reads_per_file = reads_per_file,
        read_groups = get_ReadGroups.read_groups,
        prefix = prefix,
        aligner = aligner,
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
            choices: ["mem", "aln"]
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
        disk: "10 GB"
        container: 'ghcr.io/stjudecloud/util:1.3.0'
        maxRetries: 1
    }
}
