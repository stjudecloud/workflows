version 1.1

import "../../data_structures/read_group.wdl"
import "../../tools/bwa.wdl"
import "../../tools/juicer.wdl"
import "../../tools/pairix.wdl"
import "../../tools/picard.wdl"
import "../../tools/samtools.wdl"
import "../../tools/util.wdl"
import "../general/bam-to-fastqs.wdl" as bam_to_fastqs_wf
import "../general/samtools-merge.wdl" as samtools_merge_wf

workflow hic_core {
    meta {
        description: "Hi-C core pipeline to produce unaligned BAM and .hic files"
        outputs: {
            unaligned_bam: "Unaligned BAM file",
            hic: "Juicer .hic file",
        }
        allowNestedInputs: true
    }

    parameter_meta {
        read_one_fastqs_gz: "Array of gzipped FASTQ files with 1st reads in pair. Assumes 1 file pair per read group."
        read_two_fastqs_gz: "Array of gzipped FASTQ files with 2nd reads in pair. Assumes 1 file pair per read group."
        bwa_db: "Gzipped tar archive of the bwa reference files. Files should be at the root of the archive."
        read_groups: "An array of ReadGroup objects containing the read group information to output in the BAM file. One entry per file in `read_one_fastqs_gz`/`read_two_fastqs_gz`."
        genome_id: {
            description: "Genome ID",
            choices: [
                "hg18",
                "hg19",
                "hg38",
                "dMel",
                "mm9",
                "mm10",
                "anasPlat1",
                "bTaurus3",
                "canFam3",
                "equCab2",
                "galGal4",
                "Pf3D7",
                "sacCer3",
                "sCerS288c",
                "susScr3",
                "TAIR10"
            ],
        }
        prefix: "Prefix for the BAM file. The extension `.bam` will be added."
        validate_input: "Ensure input BAM is well-formed before beginning harmonization?"
        use_all_cores: "Use all cores? Recommended for cloud environments."
        restriction_sites: {
            description: "Calculate fragment map. Requires restriction site file; each line should start with the chromosome name followed by the position of each restriction site on that chromosome, in numeric order, and ending with the size of the chromosome.",
            external_help: "https://github.com/aidenlab/juicer/wiki/Pre#restriction-site-file-format",
        }
    }

    input {
        File bwa_db
        Array[File] read_one_fastqs_gz
        Array[File] read_two_fastqs_gz
        Array[ReadGroup] read_groups
        String prefix
        String genome_id
        File? restriction_sites
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

    scatter (rg in read_groups) {
        call read_group.read_group_to_string { input: read_group = rg }
    }

    Array[String] read_groups_bwa = prefix("@RG\t", read_group_to_string.stringified_read_group)

    scatter (tuple in zip(
        zip(read_one_fastqs_gz, read_two_fastqs_gz),
        read_groups_bwa
    )) {
        call bwa.bwa_mem { input:
            read_one_fastq_gz = tuple.left.left,
            read_two_fastq_gz = tuple.left.right,
            bwa_db_tar_gz = bwa_db,
            # find tab literals, replace with '\\t' (which must be written as '\\\\t')
            # '\\t' is subbed into command blocks as '\t'
            read_group = sub(tuple.right, "\t", "\\\\t"),
            skip_mate_rescue = true,
            skip_pairing = true,
            split_smallest = true,
            short_secondary = true,
            use_all_cores,
        }
        call picard.sort { input:
            bam = bwa_mem.bam,
        }
    }
    call samtools_merge_wf.samtools_merge as merge_bam{ input:
        bams = sort.sorted_bam,
        prefix,
        use_all_cores,
    }

    call pairix.bam2pairs { input:
        bam = merge_bam.merged_bam,
        prefix,
    }

    call juicer.pre { input:
        pairs = bam2pairs.pairs,
        genome_id,
        restriction_sites,
    }

    output {
        File unaligned_bam = merge_unaligned_bam.merged_bam
        File hic = pre.hic
    }
}
