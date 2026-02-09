## **WARNING:** this workflow is experimental! Use at your own risk!
version 1.1

import "../../tools/bwa.wdl"
import "../../tools/fastp.wdl" as fp
import "../../tools/picard.wdl"
import "../../tools/samtools.wdl"
import "../../tools/util.wdl"
import "../general/samtools-merge.wdl" as samtools_merge_wf

workflow dnaseq_core_experimental {
    meta {
        name: "DNA-Seq Core (Experimental)"
        description: "Aligns DNA reads using bwa"
        warning: "We recommend against calling this workflow directly, and would suggest instead running `dnaseq_standard_experimental` or `dnaseq_standard_fastq_experimental`. Both wrapper workflows provide a nicer user experience than this workflow and will get you equivalent results."
        outputs: {
            harmonized_bam: "Harmonized DNA-Seq BAM, aligned with bwa",
            harmonized_bam_index: "Index for the harmonized DNA-Seq BAM file",
            fastp_reports: "An array of `fastp` reports (in HTML format) corresponding to each read group",
            fastp_jsons: "An array of `fastp` reports (in JSON format) corresponding to each read group",
        }
        allowNestedInputs: true
    }

    parameter_meta {
        bwa_db: "Gzipped tar archive of the bwa reference files. Files should be at the root of the archive."
        read_one_fastqs_gz: "Input gzipped FASTQ format file(s) with 1st read in pair to align"
        read_two_fastqs_gz: "Input gzipped FASTQ format file(s) with 2nd read in pair to align"
        read_groups: {
            description: "This is functionally an array of SAM `@RG` header records.",
            warning: "You should not write this input manually, but instead rely on the `ReadGroup` struct defined in `data_structures/read_group.wdl` and the utility workflow `read_group_to_string` with `format_as_sam_record = true`.",
        }
        prefix: "Prefix for the BAM file. The extension `.bam` will be added."
        aligner: {
            description: "BWA aligner to use",
            choices: [
                "mem",
                "aln",
            ],
        }
        enable_read_trimming: "Enable read trimming with `fastp`?"
        use_all_cores: "Use all cores? Recommended for cloud environments."
        reads_per_file: "Controls the number of reads per FASTQ file for internal split to run BWA in parallel."
    }

    input {
        File bwa_db
        Array[File] read_one_fastqs_gz
        Array[File] read_two_fastqs_gz
        Array[String] read_groups
        String prefix
        String aligner
        Boolean enable_read_trimming
        Boolean use_all_cores
        Int reads_per_file = 10000000
    }

    scatter (fq in read_one_fastqs_gz) {
        String read_one_names = basename(fq)
    }
    scatter (fq in read_two_fastqs_gz) {
        String read_two_names = basename(fq)
    }

    call util.check_fastq_and_rg_concordance as validate { input:
        read_one_names,
        read_two_names,
        read_groups,
    }

    scatter (tuple in zip(zip(read_one_fastqs_gz, read_two_fastqs_gz), read_groups)) {
        if (enable_read_trimming) {
            call fp.fastp as trim after validate { input:
                read_one_fastq = tuple.left.left,
                read_two_fastq = tuple.left.right,
                output_fastq = enable_read_trimming,
            }
        }
        if (!enable_read_trimming) {
            call fp.fastp after validate { input:
                read_one_fastq = tuple.left.left,
                read_two_fastq = tuple.left.right,
                output_fastq = enable_read_trimming,
            }
        }
        File chosen_r1_fastq = select_first([
            trim.read_one_fastq_gz,
            tuple.left.left,
        ])
        File chosen_r2_fastq = select_first([
            trim.read_two_fastq_gz,
            tuple.left.right,
        ])

        call util.split_fastq as read_ones after validate { input:
            fastq = chosen_r1_fastq,
            reads_per_file,
        }
        call util.split_fastq as read_twos after validate { input:
            fastq = chosen_r2_fastq,
            reads_per_file,
        }

        scatter (t in zip(read_ones.fastqs, read_twos.fastqs)) {
            if (aligner == "mem") {
                call bwa.bwa_mem { input:
                    read_one_fastq_gz = t.left,
                    read_two_fastq_gz = t.right,
                    bwa_db_tar_gz = bwa_db,
                    prefix = sub(sub(basename(t.left), "(\\.subsampled)?\\.(fastq|fq)(\\.gz)?$",
                        ""), "\\.([rR][12])\\.", "."),
                    read_group = tuple.right,
                    use_all_cores,
                }
            }
            if (aligner == "aln") {
                call bwa.bwa_aln_pe { input:
                    read_one_fastq_gz = t.left,
                    read_two_fastq_gz = t.right,
                    bwa_db_tar_gz = bwa_db,
                    prefix = sub(sub(basename(t.left), "(\\.subsampled)?\\.(fastq|fq)(\\.gz)?$",
                        ""), "\\.([rR][12])\\.", "."),
                    read_group = tuple.right,
                    use_all_cores,
                }
            }
            call picard.sort as sort { input:
                bam = select_first([
                    bwa_mem.bam,
                    bwa_aln_pe.bam,
                ]),
            }
        }
    }
    call samtools_merge_wf.samtools_merge as merge { input:
        bams = flatten(sort.sorted_bam),
        prefix,
        use_all_cores,
    }

    call samtools.index { input:
        bam = merge.merged_bam,
    }

    output {
        File harmonized_bam = merge.merged_bam
        File harmonized_bam_index = index.bam_index
        Array[File] fastp_reports = select_all(flatten([
            fastp.report,
            trim.report,
        ]))
        Array[File] fastp_jsons = select_all(flatten([
            fastp.report_json,
            trim.report_json,
        ]))
    }
}
