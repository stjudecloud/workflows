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
import "./hic-simple-core.wdl" as hic_core

workflow hic_standard {
    meta {
        description: "hi-c"
        outputs: {
            unaligned_bam: "Unaligned BAM file",
            hic: "Juicer .hic file",
        }
        allowNestedInputs: true
    }

    parameter_meta {
        bam: "Input BAM to realign"
        bwa_db: "Gzipped tar archive of the bwa reference files. Files should be at the root of the archive."
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
            help: "Common restriction sites can be downloaded from: https://bcm.app.box.com/s/19807ji76uy20cd1wau9g146nypsv2u0",
        }
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

    call parse_input { input:
        genome_id,
    }

    if (validate_input) {
        call picard.validate_bam as validate_input_bam { input:
            bam,
        }
    }

    call read_group.get_read_groups { input:
        bam,
    }

    call bam_to_fastqs_wf.bam_to_fastqs after parse_input { input:
        bam,
        paired_end = true,  # matches default but prevents user from overriding
        use_all_cores,
    }

    call hic_core.hic_core after parse_input { input:
        read_one_fastqs_gz = bam_to_fastqs.read1s,
        read_two_fastqs_gz = select_all(bam_to_fastqs.read2s),
        bwa_db,
        read_groups = get_read_groups.read_groups,
        genome_id,
        prefix,
        use_all_cores,
        restriction_sites,
    }

    output {
        File unaligned_bam = hic_core.unaligned_bam
        File hic = hic_core.hic
    }
}

task parse_input {
    meta {
        description: "Parses and validates the `hic-standard` workflow's provided inputs"
        outputs: {
            check: "Dummy output to indicate success and to enable call-caching"
        }
    }

    parameter_meta {
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
    }

    input {
        String genome_id
    }

    command <<<
        if [ $(echo ~{genome_id} | grep -Ewc "hg18|hg19|hg38|dMel|mm9|mm10|anasPlat1|bTaurus3|canFam3|equCab2|galGal4|Pf3D7|sacCer3|sCerS288c|susScr3|TAIR10") -ne 1 ]
        then
            >&2 echo "Invalid genome_id: ~{genome_id}"
            exit 1
        fi
    >>>

    output {
        String check = "passed"
    }

    runtime {
        memory: "4 GB"
        disks: "10 GB"
        container: "ghcr.io/stjudecloud/util:1.3.0"
        maxRetries: 1
    }
}
