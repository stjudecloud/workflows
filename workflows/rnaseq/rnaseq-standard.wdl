version 1.1

import "../../data_structures/read_group.wdl"
import "../../tools/picard.wdl"
import "../../tools/samtools.wdl"
import "../general/bam-to-fastqs.wdl" as bam_to_fastqs_wf
import "./rnaseq-core.wdl" as rnaseq_core_wf

workflow rnaseq_standard {
    meta {
        name: "RNA-Seq Standard"
        description: "Runs the STAR RNA-Seq alignment workflow for St. Jude Cloud"
        category: "Harmonization"
        outputs: {
            harmonized_bam: "Harmonized RNA-Seq BAM",
            bam_index: "BAI index file associated with `bam`",
            bam_checksum: "STDOUT of the `md5sum` command run on the harmonized BAM that has been redirected to a file",
            star_log: "Summary mapping statistics after mapping job is complete",
            bigwig: "BigWig format coverage file generated from `bam`",
            feature_counts: "A two column headerless TSV file. First column is feature names and second column is counts.",
            inferred_strandedness: "TSV file containing the `ngsderive strandedness` report",
            inferred_strandedness_string: "Derived strandedness from `ngsderive strandedness`",
        }
        allowNestedInputs: true
    }

    parameter_meta {
        bam: "Input BAM format file to harmonize"
        gtf: "Gzipped GTF feature file"
        star_db: "Database of reference files for the STAR aligner. The name of the root directory which was archived must match the archive's filename without the `.tar.gz` extension. Can be generated by `star-db-build.wdl`"
        contaminant_db: "A compressed reference database corresponding to the aligner chosen with `xenocp_aligner` for the contaminant genome"
        prefix: "Prefix for output files"
        xenocp_aligner: {
            description: "Aligner to use to map reads to the host genome for detecting contamination",
            choices: [
                "bwa aln",
                "bwa mem",
                "star",
            ],
        }
        strandedness: {
            description: "Strandedness protocol of the RNA-Seq experiment. If unspecified, strandedness will be inferred by `ngsderive`.",
            choices: [
                "",
                "Stranded-Reverse",
                "Stranded-Forward",
                "Unstranded",
            ],
        }
        mark_duplicates: "Add SAM flag to computationally determined duplicate reads?"
        cleanse_xenograft: "Use XenoCP to unmap reads from contaminant genome?"
        validate_input: "Ensure input BAM is well-formed before beginning harmonization?"
        use_all_cores: "Use all cores for multi-core steps?"
        subsample_n_reads: "Only process a random sampling of `n` reads. Any `n`<=`0` for processing entire input."
    }

    input {
        File bam
        File gtf
        File star_db
        File? contaminant_db
        String prefix = basename(bam, ".bam")
        String xenocp_aligner = "star"
        String strandedness = ""
        Boolean mark_duplicates = false
        Boolean cleanse_xenograft = false
        Boolean validate_input = true
        Boolean use_all_cores = false
        Int subsample_n_reads = -1
    }

    call parse_input { input:
        strand = strandedness,
        cleanse_xenograft,
        contaminant_db = defined(contaminant_db),
    }

    if (validate_input) {
        call picard.validate_bam as validate_input_bam after parse_input { input:
            bam,
        }
    }

    if (subsample_n_reads > 0) {
        call samtools.subsample after validate_input_bam { input:
            bam,
            desired_reads = subsample_n_reads,
            use_all_cores,
        }
    }
    File selected_bam = select_first([subsample.sampled_bam, bam])

    call read_group.get_read_groups after validate_input_bam { input:
        bam = selected_bam,
    }
    scatter (rg in get_read_groups.read_groups) {
        call read_group.read_group_to_string { input:
            read_group = rg,
        }
    }
    call bam_to_fastqs_wf.bam_to_fastqs after validate_input_bam { input:
        bam = selected_bam,
        paired_end = true,  # matches default but prevents user from overriding
        use_all_cores,
    }

    call rnaseq_core_wf.rnaseq_core { input:
        read_one_fastqs_gz = bam_to_fastqs.read1s,
        read_two_fastqs_gz = select_all(bam_to_fastqs.read2s),
        read_groups = read_group_to_string.validated_read_group,
        prefix,
        gtf,
        star_db,
        mark_duplicates,
        contaminant_db,
        cleanse_xenograft,
        xenocp_aligner,
        strandedness,
        use_all_cores,
    }

    output {
        File harmonized_bam = rnaseq_core.bam
        File bam_index = rnaseq_core.bam_index
        File bam_checksum = rnaseq_core.bam_checksum
        File star_log = rnaseq_core.star_log
        File bigwig = rnaseq_core.bigwig
        File feature_counts = rnaseq_core.feature_counts
        File inferred_strandedness = rnaseq_core.inferred_strandedness
        String inferred_strandedness_string = rnaseq_core.inferred_strandedness_string
    }
}

task parse_input {
    meta {
        description: "Parses and validates the `rnaseq_standard[_fastq]` workflows' provided inputs"
    }

    parameter_meta {
        strand: {
            description: "Provided strandedness protocol of the RNA-Seq experiment",
            choices: [
                "",
                "Stranded-Reverse",
                "Stranded-Forward",
                "Unstranded",
            ],
        }
        cleanse_xenograft: "Use XenoCP to unmap reads from contaminant genome?"
        contaminant_db: "Was a `contaminant_db` supplied by the user? Must `true` if `cleanse_xenograft` is `true`."
    }

    input {
        String strand
        Boolean cleanse_xenograft
        Boolean contaminant_db
    }

    command <<<
        if [ -n "~{strand}" ] \
            && [ "~{strand}" != "Stranded-Reverse" ] \
            && [ "~{strand}" != "Stranded-Forward" ] \
            && [ "~{strand}" != "Unstranded" ]
        then
            >&2 echo "strandedness must be:"
            >&2 echo "'', 'Stranded-Reverse', 'Stranded-Forward', or 'Unstranded'"
            exit 1
        fi

        if ~{cleanse_xenograft} && ! ~{contaminant_db}; then
            >&2 echo "'contaminant_db' must be supplied if 'cleanse_xenograft' is 'true'"
            exit 1
        fi
    >>>

    runtime {
        container: "ghcr.io/stjudecloud/util:2.2.2"
        maxRetries: 1
    }
}
