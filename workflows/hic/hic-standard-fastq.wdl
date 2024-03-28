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
        read_one_fastqs_gz: "Array of gzipped FASTQ files with 1st reads in pair"
        read_two_fastqs_gz: "Array of gzipped FASTQ files with 2nd reads in pair"
        bwa_db: "Gzipped tar archive of the bwa reference files. Files should be at the root of the archive."
        read_groups: {
            description: "An Array of structs defining read groups to include in the harmonized BAM. Must correspond to input FASTQs. Each read group ID must be contained in the basename of a FASTQ file or pair of FASTQ files if Paired-End. This requirement means the length of `read_groups` must equal the length of `read_one_fastqs_gz` and the length of `read_two_fastqs_gz` if non-zero. Only the `ID` field is required, and it must be unique for each read group defined. See top of file for help formatting your input JSON.",  # TODO handle unknown RG case
            external_help: "https://samtools.github.io/hts-specs/SAMv1.pdf",
            fields: {
                ID: "Read group identifier. Each Read Group must have a unique ID. The value of ID is used in the RG tags of alignment records.",
                BC: "Barcode sequence identifying the sample or library. This value is the expected barcode bases as read by the sequencing machine in the absence of errors. If there are several barcodes for the sample/library (e.g., one on each end of the template), the recommended implementation concatenates all the barcodes separating them with hyphens (`-`).",
                CN: "Name of sequencing center producing the read.",
                DS: "Description.",
                DT: "Date the run was produced (ISO8601 date or date/time).",
                FO: "Flow order. The array of nucleotide bases that correspond to the nucleotides used for each flow of each read. Multi-base flows are encoded in IUPAC format, and non-nucleotide flows by various other characters. Format: /\\*|[ACMGRSVTWYHKDBN]+/",
                KS: "The array of nucleotide bases that correspond to the key sequence of each read.",
                LB: "Library.",
                PG: "Programs used for processing the read group.",
                PI: "Predicted median insert size, rounded to the nearest integer.",
                PL: "Platform/technology used to produce the reads. Valid values: CAPILLARY, DNBSEQ (MGI/BGI), ELEMENT, HELICOS, ILLUMINA, IONTORRENT, LS454, ONT (Oxford Nanopore), PACBIO (Pacific Biosciences), SINGULAR, SOLID, and ULTIMA. This field should be omitted when the technology is not in this list (though the PM field may still be present in this case) or is unknown.",
                PM: "Platform model. Free-form text providing further details of the platform/technology used.",
                PU: "Platform unit (e.g., flowcell-barcode.lane for Illumina or slide for SOLiD). Unique identifier.",
                SM: "Sample. Use pool name where a pool is being sequenced."
            }
        }
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
        String genomeID = "hg38"
        String prefix = basename(read_one_fastqs_gz[0], ".fq.gz")
        Boolean validate_input = true
        Boolean use_all_cores = false
    }

    call hic_core.hic_core { input:
        read_one_fastqs_gz=read_one_fastqs_gz,
        read_two_fastqs_gz=read_two_fastqs_gz,
        bwa_db=bwa_db,
        read_groups=read_groups,
        genomeID=genomeID,
        prefix=prefix,
        use_all_cores=use_all_cores,
    }   

    output {
        File unaligned_bam = hic_core.unaligned_bam
        File hic = hic_core.hic
    }
}