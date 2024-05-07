## **WARNING:** this workflow is experimental! Use at your own risk!

version 1.1

import "../../tools/bwa.wdl"
import "../../tools/picard.wdl"
import "../../tools/samtools.wdl"
import "../../tools/util.wdl"
import "../general/bam-to-fastqs.wdl" as bam_to_fastqs_wf
import "../general/samtools-merge.wdl" as samtools_merge_wf

workflow dnaseq_standard_experimental {
    meta {
        description: "Aligns DNA reads using bwa"
        outputs: {
            harmonized_bam: "Harmonized DNA-Seq BAM, aligned with bwa"
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

    call util.get_read_groups { input:
        bam=bam,
        format_for_star=false,
    }  # TODO what happens if no RG records?
    call bam_to_fastqs_wf.bam_to_fastqs { input:
        bam=bam,
        paired_end=true,  # matches default but prevents user from overriding
        use_all_cores=use_all_cores,
    }

    # For each read group
    scatter (tuple in zip(
        zip(bam_to_fastqs.read1s, bam_to_fastqs.read2s),
        get_read_groups.read_groups
    )) {
        call util.split_fastq as read_ones { input:
            fastq = tuple.left.left,
            reads_per_file = reads_per_file
        }

        call util.split_fastq as read_twos { input:
            fastq = select_first([tuple.left.right, "undefined"]),
            reads_per_file = reads_per_file
        }

        scatter (t in zip(read_ones.fastqs, read_twos.fastqs)) {
            if (aligner == "mem") {
                call bwa.bwa_mem { input:
                    read_one_fastq_gz=t.left,
                    read_two_fastq_gz=select_first([t.right, "undefined"]),
                    bwa_db_tar_gz=bwa_db,
                    # find tab literals, replace with '\\t' (which must be written as '\\\\t')
                    # '\\t' is subbed into command blocks as '\t'
                    read_group=sub(tuple.right, "\t", "\\\\t"),
                    use_all_cores=use_all_cores,
                }
            }
            if ( aligner == "aln") {
                call bwa.bwa_aln_pe { input:
                    read_one_fastq_gz=t.left,
                    read_two_fastq_gz=select_first([t.right, "undefined"]),
                    bwa_db_tar_gz=bwa_db,
                    # find tab literals, replace with '\\t' (which must be written as '\\\\t')
                    # '\\t' is subbed into command blocks as '\t'
                    read_group=sub(tuple.right, "\t", "\\\\t"),
                    use_all_cores=use_all_cores,
                }
            }
            call picard.sort as sort { input:
                 bam=select_first([bwa_mem.bam, bwa_aln_pe.bam])
            }

        }

        call samtools_merge_wf.samtools_merge as inner_merge { input:
            bams = sort.sorted_bam,
            prefix = basename(tuple.left.left, ".fastq.gz"),
            use_all_cores = use_all_cores,
        }
    }
    call samtools_merge_wf.samtools_merge as rg_merge { input:
        bams=inner_merge.merged_bam,
        prefix=prefix,
        use_all_cores=use_all_cores,
    }

    call samtools.index { input:
        bam=rg_merge.merged_bam,
    }

    output {
        File harmonized_bam = rg_merge.merged_bam
        File harmonized_bam_index = index.bam_index
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
