## **WARNING:** this workflow is experimental! Use at your own risk!

version 1.1

import "../../data_structures/read_group.wdl"
import "../../tools/bwa.wdl"
import "../../tools/picard.wdl"
import "../../tools/samtools.wdl"
import "../../tools/util.wdl"
import "../general/samtools-merge.wdl" as samtools_merge_wf

workflow dnaseq_core_experimental {
    meta {
        name: "DNA-Seq Core (Experimental)"
        description: "Aligns DNA reads using bwa"
        outputs: {
            harmonized_bam: "Harmonized DNA-Seq BAM, aligned with bwa",
            harmonized_bam_index: "Index for the harmonized DNA-Seq BAM file",
        }
        allowNestedInputs: true
    }

    parameter_meta {
        bwa_db: "Gzipped tar archive of the bwa reference files. Files should be at the root of the archive."
        read_one_fastqs_gz: "Input gzipped FASTQ format file(s) with 1st read in pair to align"
        read_two_fastqs_gz: "Input gzipped FASTQ format file(s) with 2nd read in pair to align"
        read_groups: "An Array of structs defining read groups to include in the harmonized BAM. Must correspond to input FASTQs. Each read group ID must be contained in the basename of a pair of FASTQ files. This requirement means the length of `read_groups` must equal the length of `read_one_fastqs_gz` and the length of `read_two_fastqs_gz`. Only the `ID` field is required, and it must be unique for each read group defined. See data_structures/read_group.wdl for help formatting your input JSON."
        prefix: "Prefix for the BAM file. The extension `.bam` will be added."
        sample_override: "Value to override the SM field of *every* read group."
        aligner: {
            description: "BWA aligner to use",
            choices: [
                "mem",
                "aln"
            ],
        }
        use_all_cores: "Use all cores? Recommended for cloud environments."
        reads_per_file: "Controls the number of reads per FASTQ file for internal split to run BWA in parallel."
    }

    input {
        File bwa_db
        Array[File] read_one_fastqs_gz
        Array[File] read_two_fastqs_gz
        Array[ReadGroup] read_groups
        String prefix
        String? sample_override
        String aligner = "mem"
        Boolean use_all_cores = false
        Int reads_per_file = 10000000
    }

    scatter (tuple in zip(
        zip(read_one_fastqs_gz, read_two_fastqs_gz),
        read_groups
    )) {
        if (defined(sample_override)) {
            # override the SM field of every read group
            ReadGroup rg = ReadGroup{
                ID: tuple.right.ID,
                BC: tuple.right.BC,
                CN: tuple.right.CN,
                DS: tuple.right.DS,
                DT: tuple.right.DT,
                FO: tuple.right.FO,
                KS: tuple.right.KS,
                LB: tuple.right.LB,
                PG: tuple.right.PG,
                PI: tuple.right.PI,
                PL: tuple.right.PL,
                PM: tuple.right.PM,
                PU: tuple.right.PU,
                SM: sample_override,
            }
        }

        call read_group.read_group_to_string { input:
            read_group = select_first([rg, tuple.right])
        }

        String rg_string = "@RG " + read_group_to_string.stringified_read_group

        call util.split_fastq as read_ones { input:
            fastq = tuple.left.left,
            reads_per_file,
        }

        call util.split_fastq as read_twos { input:
            fastq = tuple.left.right,
            reads_per_file,
        }

        scatter (t in zip(read_ones.fastqs, read_twos.fastqs)) {
            if (aligner == "mem") {
                call bwa.bwa_mem { input:
                    read_one_fastq_gz = t.left,
                    read_two_fastq_gz = t.right,
                    bwa_db_tar_gz = bwa_db,
                    prefix = sub(sub(
                        basename(t.left),
                        "(\\.subsampled)?\\.(fastq|fq)(\\.gz)?$",
                        ""
                    ), "\\.([rR][12])\\.", "."),
                    # find spaces, replace with '\\t'
                    # (which must be written as '\\\\t')
                    # '\\t' is subbed into command blocks as '\t'
                    read_group = sub(rg_string, " ", "\\\\t"),
                    use_all_cores,
                }
            }
            if (aligner == "aln") {
                call bwa.bwa_aln_pe { input:
                    read_one_fastq_gz = t.left,
                    read_two_fastq_gz = t.right,
                    bwa_db_tar_gz = bwa_db,
                    prefix = sub(sub(
                        basename(t.left),
                        "(\\.subsampled)?\\.(fastq|fq)(\\.gz)?$",
                        ""
                    ), "\\.([rR][12])\\.", "."),
                    # find spaces, replace with '\\t'
                    # (which must be written as '\\\\t')
                    # '\\t' is subbed into command blocks as '\t'
                    read_group = sub(rg_string, " ", "\\\\t"),
                    use_all_cores,
                }
            }
            call picard.sort as sort { input:
                 bam = select_first([bwa_mem.bam, bwa_aln_pe.bam])
            }
        }
    }
    call samtools_merge_wf.samtools_merge as rg_merge { input:
        bams = flatten(sort.sorted_bam),
        prefix,
        use_all_cores,
    }

    call samtools.index { input:
        bam = rg_merge.merged_bam,
    }

    output {
        File harmonized_bam = rg_merge.merged_bam
        File harmonized_bam_index = index.bam_index
    }
}
