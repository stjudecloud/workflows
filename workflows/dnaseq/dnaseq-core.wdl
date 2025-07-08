## **WARNING:** this workflow is experimental! Use at your own risk!

version 1.1

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
        read_groups: "TODO"
        prefix: "Prefix for the BAM file. The extension `.bam` will be added."
        aligner: {
            description: "BWA aligner to use",
            choices: [
                "mem",
                "aln",
            ],
        }
        use_all_cores: "Use all cores? Recommended for cloud environments."
        reads_per_file: "Controls the number of reads per FASTQ file for internal split to run BWA in parallel."
    }

    input {
        File bwa_db
        Array[File] read_one_fastqs_gz
        Array[File] read_two_fastqs_gz
        Array[String] read_groups
        String prefix
        String aligner = "mem"
        Boolean use_all_cores = false
        Int reads_per_file = 10000000
    }

    scatter (fq in read_one_fastqs_gz) {
        String read_one_basenames = basename(fq)
    }
    scatter (fq in read_two_fastqs_gz) {
        String read_two_basenames = basename(fq)
    }

    #@ except: UnusedCall
    call util.check_fastq_and_rg_concordance { input:
        read_one_basenames,
        read_two_basenames,
        read_groups,
    }

    scatter (tuple in zip(
        zip(read_one_fastqs_gz, read_two_fastqs_gz),
        read_groups
    )) {
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
                    read_group = tuple.right,
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
                    read_group = tuple.right,
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
