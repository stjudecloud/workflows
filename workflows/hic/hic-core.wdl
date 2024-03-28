version 1.1

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
            unaligned_bam: "Unaligned BAM file"
            hic: "Juicer .hic file"
        }
        allowNestedInputs: true
    }

    parameter_meta {
        read_one_fastqs_gz: "Array of gzipped FASTQ files with 1st reads in pair"
        read_two_fastqs_gz: "Array of gzipped FASTQ files with 2nd reads in pair"
        bwa_db: "Gzipped tar archive of the bwa reference files. Files should be at the root of the archive."
        read_groups: "A string containing the read group information to output in the BAM file. If including multiple read group fields per-read group, they should be space delimited. Read groups should be comma separated, with a space on each side (i.e. ' , '). The ID field must come first for each read group and must be contained in the basename of a FASTQ file or pair of FASTQ files if Paired-End. Example: `ID:rg1 PU:flowcell1.lane1 SM:sample1 PL:illumina LB:sample1_lib1 , ID:rg2 PU:flowcell1.lane2 SM:sample1 PL:illumina LB:sample1_lib1`. These two read groups could be associated with the following four FASTQs: `sample1.rg1_R1.fastq,sample1.rg2_R1.fastq` and `sample1.rg1_R2.fastq,sample1.rg2_R2.fastq`"
        genomeID: {
            description: "Genome ID"
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
            ]
        }
        prefix: "Prefix for the BAM file. The extension `.bam` will be added."
        validate_input: "Ensure input BAM is well-formed before beginning harmonization?"
        use_all_cores: "Use all cores? Recommended for cloud environments."
    }

    input {
        Array[File] read_one_fastqs_gz
        Array[File] read_two_fastqs_gz
        File bwa_db
        Array[ReadGroup] read_groups
        String prefix
        String genomeID = "hg38"
        Boolean validate_input = true
        Boolean use_all_cores = false
    }

    scatter (tuple in zip(
        zip(read_one_fastqs_gz, read_two_fastqs_gz),
        read_groups
    )) {
        call picard.fastq_to_sam {
            input:
                read_one_fastq_gz=tuple.left.left,
                read_two_fastq_gz=tuple.left.right,
                read_group=tuple.right
                prefix=prefix,
                read_group_name=tuple.right.ID,
                sample_name=tuple.right.SM,
                library_name=tuple.right.LB,
                sequencing_center=tuple.right.CN,
                run_date=tuple.right.DT,
                platform_unit=tuple.right.PU,
                platform=tuple.right.PL,
                platform_model=tuple.right.PM,
        }
    }

    call samtools_merge_wf.samtools_merge as merge_unaligned_bam { input:
        bams=fastq_to_sam.unaligned_bam,
        prefix=prefix,
        use_all_cores=use_all_cores,
    }

    scatter (tuple in zip(
        zip(read_one_fastqs_gz, read_two_fastqs_gz),
        read_groups
    )) {
        call bwa.bwa_mem { input:
            read_one_fastq_gz=tuple.left.left,
            read_two_fastq_gz=tuple.left.right,
            bwa_db_tar_gz=bwa_db,
            # find tab literals, replace with '\\t' (which must be written as '\\\\t')
            # '\\t' is subbed into command blocks as '\t'
            read_group=sub(tuple.right, "\t", "\\\\t"),
            skip_mate_rescue=true,
            skip_pairing=true,
            split_smallest=true,
            short_secondary=true,
            use_all_cores=use_all_cores,
        }
        call picard.sort { input:
            bam=bwa_mem.bam,
        }
    }
    call samtools_merge_wf.samtools_merge as merge_bam{ input:
        bams=sort.sorted_bam,
        prefix=prefix,
        use_all_cores=use_all_cores,
    }
    
    call pairix.bam2pairs {
        input:
            bam=merge_bam.merged_bam,
            prefix=prefix,
    }

    call juicer.pre {
        input:
            pairs=bam2pairs.pairs,
            genomeID=genomeID,
    }

    output {
        File unaligned_bam = "x.bam" 
        File hic = pre.hic
    }
}