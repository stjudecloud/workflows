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
import "./hic-core.wdl" as hic_core

workflow hic_standard {
    meta {
        description: "hi-c"
        outputs: {
            unaligned_bam: "Unaligned BAM file"
            hic: "Juicer .hic file"
        }
        allowNestedInputs: true
    }
    
    parameter_meta {
        bam: "Input BAM to realign"
        bwa_db: "Gzipped tar archive of the bwa reference files. Files should be at the root of the archive."
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
        File bam
        File bwa_db
        String genomeID = "hg38"
        String prefix = basename(bam, ".bam")
        Boolean validate_input = true
        Boolean use_all_cores = false
    }

    if (validate_input) {
        call picard.validate_bam as validate_input_bam { input:
            bam=bam,
        }
    }

    call read_group.get_ReadGroups { input:
        bam=bam,
    }

    call bam_to_fastqs_wf.bam_to_fastqs { input:
        bam=bam,
        paired_end=true,  # matches default but prevents user from overriding
        use_all_cores=use_all_cores,
    }

    call hic_core.hic_core { input:
        read_one_fastqs_gz=bam_to_fastqs.read1s,
        read_two_fastqs_gz=select_all(bam_to_fastqs.read2s),
        bwa_db=bwa_db,
        read_groups=get_ReadGroups.read_groups,
        genomeID=genomeID,
        prefix=prefix,
        use_all_cores=use_all_cores,
    }   

    output {
        File unaligned_bam = hic_core.unaligned_bam
        File hic = hic_core.hic
    }
}