# SPDX-License-Identifier: MIT
# Copyright St. Jude Children's Research Hospital
version 1.1

import "../../tools/fastqc.wdl" as fastqc_tasks
import "../../tools/fq.wdl"
import "../../tools/kraken2.wdl"
import "../../tools/md5sum.wdl"
import "../../tools/mosdepth.wdl"
import "../../tools/multiqc.wdl" as multiqc_tasks
import "../../tools/ngsderive.wdl"
import "../../tools/picard.wdl"
import "../../tools/qualimap.wdl"
import "../../tools/samtools.wdl"
import "../../tools/util.wdl"
import "./markdups-post.wdl" as markdups_post_wf

workflow quality_check {
    meta {
        description: "Performs comprehensive quality checks, aggregating all analyses and metrics into a final MultiQC report."
        external_help: "https://multiqc.info/"
        outputs: {
            bam_checksum: "STDOUT of the `md5sum` command run on the input BAM that has been redirected to a file"
            validate_sam_file: "Validation report produced by `picard ValidateSamFile`. Validation warnings and errors are logged."
            mark_duplicates_metrics: {
                description: "The METRICS_FILE result of `picard MarkDuplicates`"
                external_help: "http://broadinstitute.github.io/picard/picard-metric-definitions.html#DuplicationMetrics"
            }
            flagstat_report: "`samtools flagstat` STDOUT redirected to a file. If `mark_duplicates` is `true`, then this result will be generated from the duplicate marked BAM."
            fastqc_results: "A gzipped tar archive of all FastQC output files"
            instrument_file: "TSV file containing the `ngsderive isntrument` report for the input BAM file"
            read_length_file: "TSV file containing the `ngsderive readlen` report for the input BAM file"
            inferred_encoding: "TSV file containing the `ngsderive encoding` report for the input BAM file"
            alignment_metrics: {
                description: "The text file output of `picard CollectAlignmentSummaryMetrics`"
                external_help: "http://broadinstitute.github.io/picard/picard-metric-definitions.html#AlignmentSummaryMetrics"
            }
            alignment_metrics_pdf: "The PDF file output of `picard CollectAlignmentSummaryMetrics`"
            insert_size_metrics: {
                description: "The text file output of `picard CollectInsertSizeMetrics`. If `mark_duplicates` is `true`, then this result will be generated from the duplicate marked BAM."
                external_help: "http://broadinstitute.github.io/picard/picard-metric-definitions.html#InsertSizeMetrics"
            }
            insert_size_metrics_pdf: "The PDF file output of `picard CollectInsertSizeMetrics`. If `mark_duplicates` is `true`, then this result will be generated from the duplicate marked BAM."
            quality_score_distribution_txt: "The text file output of `picard QualityScoreDistribution`"
            quality_score_distribution_pdf: "The PDF file output of `picard QualityScoreDistribution`"
            phred_scores: "Headered TSV file containing PHRED score statistics"
            kraken_report: {
                description: "A Kraken2 summary report"
                external_help: "https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#sample-report-output-format"
            }
            mosdepth_global_dist: "The `$prefix.mosdepth.global.dist.txt` file contains a cumulative distribution indicating the proportion of total bases that were covered for at least a given coverage value. It does this for each chromosome, and for the whole genome."
            mosdepth_global_summary: "A summary of mean depths per chromosome"
            mosdepth_region_dist: "The `$prefix.mosdepth.region.dist.txt` file contains a cumulative distribution indicating the proportion of total bases in the region(s) defined by the `coverage_bed` that were covered for at least a given coverage value. There will be one file in this array for each `coverage_beds` input file."
            mosdepth_region_summary: "A summary of mean depths per chromosome and within specified regions per chromosome. There will be one file in this array for each `coverage_beds` input file."
            multiqc_report: "A gzipped tar archive of all MultiQC output files"
            orig_read_count: "A TSV report containing the original read count before subsampling"
            kraken_sequences: {
                description: "Detailed Kraken2 output that has been gzipped"
                external_help: "https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#standard-kraken-output-format"
            }
            mosdepth_dups_marked_global_dist: "The `$prefix.mosdepth.global.dist.txt` file contains a cumulative distribution indicating the proportion of total bases that were covered for at least a given coverage value. It does this for each chromosome, and for the whole genome. This file is produced from analyzing the duplicate marked BAM (only present if `mark_duplicates = true`)."
            mosdepth_dups_marked_global_summary: "A summary of mean depths per chromosome. This file is produced from analyzing the duplicate marked BAM (only present if `mark_duplicates = true`)."
            mosdepth_dups_marked_region_dist: "The `$prefix.mosdepth.region.dist.txt` file contains a cumulative distribution indicating the proportion of total bases in the region(s) defined by the `coverage_bed` that were covered for at least a given coverage value. There will be one file in this array for each `coverage_beds` input file. This file is produced from analyzing the duplicate marked BAM (only present if `mark_duplicates = true`)."
            mosdepth_dups_marked_region_summary: "A summary of mean depths per chromosome and within specified regions per chromosome. There will be one file in this array for each `coverage_beds` input file. This file is produced from analyzing the duplicate marked BAM (only present if `mark_duplicates = true`)."
            inferred_strandedness: "TSV file containing the `ngsderive strandedness` report"
            qualimap_rnaseq_results: "Gzipped tar archive of all QualiMap output files"
            junction_summary: "TSV file containing the `ngsderive junction-annotation` summary"
            junctions: "TSV file containing a detailed list of annotated junctions"
            IntermediateFiles: "Any and all files produced as intermediate during pipeline processing. Only output if `output_intermediate_files = true`."
        }
    }

    parameter_meta {
        bam: "Input BAM format file to quality check"
        bam_index: "BAM index file corresponding to the input BAM"
        kraken_db: "Kraken2 database. Can be generated with `make-qc-reference.wdl`. Must be a tarball without a root directory."
        molecule: {
            description: "Data type"
            choices: [
                'DNA',
                'RNA'
            ]
        }
        gtf: "GTF features file. Gzipped or uncompressed. **Required** for RNA-Seq data."
        multiqc_config: "YAML file for configuring MultiQC"
        extra_multiqc_inputs: "An array of additional files to pass directly into MultiQC"
        coverage_beds: "An array of 3 column BEDs which are passed to the `-b` flag of mosdepth, in order to restrict coverage analysis to select regions"
        coverage_labels: "An array of equal length to `coverage_beds` which determines the prefix label applied to the output files. If omitted, defaults of `regions1`, `regions2`, etc. will be used."
        prefix: "Prefix for all results files"
        mark_duplicates: "Mark duplicates before analyses? Note that regardless of this setting, `picard MarkDuplicates` will be run in order to generate a `*.MarkDuplicates.metrics.txt` file. However if `mark_duplicates` is set to `false`, no BAM will be generated. If set to `true`, a BAM will be generated and passed to selected downstream analyses."
        output_intermediate_files: "Output intermediate files? FASTQs, if RNA a collated BAM, if `mark_duplicates==true` a duplicate marked BAM with an index and MD5. *WARNING* these files can be large."
        use_all_cores: "Use all cores? Recommended for cloud environments. Not recommended for cluster environments."
        subsample_n_reads: "Only process a random sampling of `n` reads. Any `n`<=`0` for processing entire input. Subsampling is done probabalistically so the exact number of reads in the output will have some variation."
    }

    input {
        File bam
        File bam_index
        File kraken_db
        String molecule
        File? gtf
        File multiqc_config = "https://raw.githubusercontent.com/stjudecloud/workflows/main/workflows/qc/inputs/multiqc_config_hg38.yaml"
        Array[File] extra_multiqc_inputs = []
        Array[File] coverage_beds = []
        Array[String] coverage_labels = []
        String prefix = basename(bam, ".bam")
        Boolean mark_duplicates = molecule == "RNA"
        Boolean output_intermediate_files = false
        Boolean use_all_cores = false
        Int subsample_n_reads = -1
    }

    call parse_input { input:
        gtf_provided=defined(gtf),
        input_molecule=molecule,
        coverage_beds_len=length(coverage_beds),
        coverage_labels=coverage_labels,
    }

    call md5sum.compute_checksum { input: file=bam }

    call samtools.quickcheck { input: bam=bam }
    call util.compression_integrity { input: bgzipped_file=bam }

    if (subsample_n_reads > 0) {
        call samtools.subsample { input:
            bam=quickcheck.checked_bam,
            prefix=prefix,
            desired_reads=subsample_n_reads,
            use_all_cores=use_all_cores,
        }
        if (defined(subsample.sampled_bam)) {
            call samtools.index as subsample_index { input:
                bam=select_first([subsample.sampled_bam, "undefined"]),
                use_all_cores=use_all_cores,
            }
        }
    }
    # If subsampling is disabled or input BAM has fewer reads than `subsample_n_reads`
    # this will be `quickcheck.checked_bam`
    File post_subsample_bam = select_first([
        subsample.sampled_bam,
        quickcheck.checked_bam
    ])
    File post_subsample_bam_index = select_first([
        subsample_index.bam_index,
        bam_index
    ])
    String post_subsample_prefix = if (defined(subsample.sampled_bam))
        then prefix + ".subsampled"
        else prefix

    call picard.validate_bam { input:
        bam=post_subsample_bam,
        outfile_name=post_subsample_prefix + ".ValidateSamFile.txt",
        succeed_on_errors=true,
        ignore_list=[],
        summary_mode=true,
    }

    call picard.collect_alignment_summary_metrics { input:
        bam=post_subsample_bam,
        prefix=post_subsample_prefix + ".CollectAlignmentSummaryMetrics",
    }
    call picard.quality_score_distribution { input:
        bam=post_subsample_bam,
        prefix=post_subsample_prefix + ".QualityScoreDistribution",
    }
    call fastqc_tasks.fastqc { input:
        bam=post_subsample_bam,
        prefix=post_subsample_prefix + ".fastqc_results",
        use_all_cores=use_all_cores,
    }
    call ngsderive.instrument { input:
        bam=post_subsample_bam,
        outfile_name=post_subsample_prefix + ".instrument.tsv",
    }
    call ngsderive.read_length { input:
        bam=post_subsample_bam,
        bam_index=post_subsample_bam_index,
        outfile_name=post_subsample_prefix + ".readlength.tsv",
    }
    call ngsderive.encoding { input:
        ngs_files=[post_subsample_bam],
        outfile_name=post_subsample_prefix + ".encoding.tsv",
        num_reads=-1,
    }
    call ngsderive.endedness { input:
        bam=post_subsample_bam,
        outfile_name=post_subsample_prefix + ".endedness.tsv",
        lenient=true,
    }
    call util.global_phred_scores { input:
        bam=post_subsample_bam,
        prefix=post_subsample_prefix,
    }

    call samtools.collate_to_fastq { input:
        bam=post_subsample_bam,
        prefix=post_subsample_prefix,
        # RNA needs a collated BAM for Qualimap
        # DNA can skip the associated storage costs
        store_collated_bam=(molecule == "RNA"),
        # disabling fast_mode enables writing of secondary and supplementary alignments
        # to the collated BAM when processing RNA.
        # Those alignments are used downstream by Qualimap.
        fast_mode=(molecule != "RNA"),
        paired_end=true,  # matches default but prevents user from overriding
        interleaved=false,  # matches default but prevents user from overriding
        use_all_cores=use_all_cores,
    }

    call fq.fqlint { input:
        read_one_fastq=select_first([collate_to_fastq.read_one_fastq_gz, "undefined"]),
        read_two_fastq=collate_to_fastq.read_two_fastq_gz,
    }
    call kraken2.kraken { input:
        read_one_fastq_gz=fqlint.validated_read1,
        read_two_fastq_gz=select_first([fqlint.validated_read2, "undefined"]),
        db=kraken_db,
        prefix=post_subsample_prefix,
        use_all_cores=use_all_cores,
    }

    call mosdepth.coverage as wg_coverage { input:
        bam=post_subsample_bam,
        bam_index=post_subsample_bam_index,
        prefix=post_subsample_prefix + ".whole_genome",
    }
    scatter(coverage_pair in zip(coverage_beds, parse_input.labels)) {
        call mosdepth.coverage as regions_coverage { input:
            bam=post_subsample_bam,
            bam_index=post_subsample_bam_index,
            coverage_bed=coverage_pair.left,
            prefix=post_subsample_prefix + "." + coverage_pair.right,
        }
    }

    if (molecule == "RNA") {
        call ngsderive.junction_annotation { input:
            bam=post_subsample_bam,
            bam_index=post_subsample_bam_index,
            gene_model=select_first([gtf, "undefined"]),
            prefix=post_subsample_prefix,
        }
        call ngsderive.strandedness { input:
            bam=post_subsample_bam,
            bam_index=post_subsample_bam_index,
            gene_model=select_first([gtf, "undefined"]),
            outfile_name=post_subsample_prefix + ".strandedness.tsv",
        }
        call qualimap.rnaseq as qualimap_rnaseq { input:
            bam=select_first([collate_to_fastq.collated_bam, "undefined"]),
            prefix=post_subsample_prefix + ".qualimap_rnaseq_results",
            gtf=select_first([gtf, "undefined"]),
            name_sorted=true,
            paired_end=true,  # matches default but prevents user from overriding
        }
    }

    call picard.mark_duplicates as markdups { input:
        bam=post_subsample_bam,
        create_bam=mark_duplicates,
        prefix=post_subsample_prefix + ".MarkDuplicates",
    }
    if (mark_duplicates) {
        call markdups_post_wf.markdups_post { input:
            markdups_bam=select_first([
                markdups.duplicate_marked_bam,
                "undefined"
            ]),
            markdups_bam_index=select_first([
                markdups.duplicate_marked_bam_index,
                "undefined"
            ]),
            coverage_beds=coverage_beds,
            coverage_labels=parse_input.labels,
            prefix=post_subsample_prefix + ".MarkDuplicates",
        }
    }
    if (! mark_duplicates) {
        # These analyses are called in the markdups_post workflow.
        # They should still be run if duplicates were not marked.
        call picard.collect_insert_size_metrics { input:
            bam=post_subsample_bam,
            prefix=post_subsample_prefix + ".CollectInsertSizeMetrics",
        }
        call samtools.flagstat { input:
            bam=post_subsample_bam,
            outfile_name=post_subsample_prefix + ".flagstat.txt",
        }
    }

    call multiqc_tasks.multiqc { input:
        input_files=select_all(flatten([
            [
                validate_bam.validate_report,
                markdups.mark_duplicates_metrics,
                flagstat.flagstat_report,
                markdups_post.flagstat_report,
                instrument.instrument_file,
                read_length.read_length_file,
                encoding.encoding_file,
                endedness.endedness_file,
                fastqc.raw_data,
                collect_alignment_summary_metrics.alignment_metrics,
                collect_insert_size_metrics.insert_size_metrics,
                markdups_post.insert_size_metrics,
                quality_score_distribution.quality_score_distribution_txt,
                kraken.report,
                wg_coverage.summary,
                wg_coverage.global_dist,
                global_phred_scores.phred_scores,
                subsample.orig_read_count,
                markdups_post.mosdepth_global_summary,
                markdups_post.mosdepth_global_dist,
                strandedness.strandedness_file,
                junction_annotation.junction_summary,
                qualimap_rnaseq.raw_summary,
                qualimap_rnaseq.raw_coverage
            ],
            regions_coverage.summary,
            regions_coverage.region_dist,
            extra_multiqc_inputs,
            select_first([markdups_post.mosdepth_region_summary, []]),
            select_first([markdups_post.mosdepth_region_dist, []])
        ])),
        config=multiqc_config,
        prefix=post_subsample_prefix + ".multiqc",
    }

    if (output_intermediate_files) {
        IntermediateFiles optional_files = {
            "subsampled_bam": subsample.sampled_bam,
            "subsampled_bam_index": subsample_index.bam_index,
            "collated_bam": collate_to_fastq.collated_bam,
            "read_one_fastq_gz": collate_to_fastq.read_one_fastq_gz,
            "read_two_fastq_gz": collate_to_fastq.read_two_fastq_gz,
            "singleton_reads_fastq_gz": collate_to_fastq.singleton_reads_fastq_gz,
            "interleaved_reads_fastq_gz": collate_to_fastq.interleaved_reads_fastq_gz,
            "duplicate_marked_bam": markdups.duplicate_marked_bam,
            "duplicate_marked_bam_index": markdups.duplicate_marked_bam_index,
            "duplicate_marked_bam_md5": markdups.duplicate_marked_bam_md5
        }
    }

    output {
        File bam_checksum = compute_checksum.md5sum
        File validate_sam_file = validate_bam.validate_report
        File mark_duplicates_metrics = markdups.mark_duplicates_metrics
        File flagstat_report = select_first([
            markdups_post.flagstat_report,
            flagstat.flagstat_report
        ])
        File fastqc_results = fastqc.results
        File instrument_file = instrument.instrument_file
        File read_length_file = read_length.read_length_file
        File inferred_encoding = encoding.encoding_file
        File inferred_endedness = endedness.endedness_file
        File alignment_metrics = collect_alignment_summary_metrics.alignment_metrics
        File alignment_metrics_pdf
            = collect_alignment_summary_metrics.alignment_metrics_pdf
        File insert_size_metrics = select_first([
            markdups_post.insert_size_metrics,
            collect_insert_size_metrics.insert_size_metrics
        ])
        File insert_size_metrics_pdf = select_first([
            markdups_post.insert_size_metrics_pdf,
            collect_insert_size_metrics.insert_size_metrics_pdf
        ])
        File quality_score_distribution_txt
            = quality_score_distribution.quality_score_distribution_txt
        File quality_score_distribution_pdf
            = quality_score_distribution.quality_score_distribution_pdf
        File phred_scores = global_phred_scores.phred_scores
        File kraken_report = kraken.report
        File mosdepth_global_dist = wg_coverage.global_dist
        File mosdepth_global_summary = wg_coverage.summary
        Array[File] mosdepth_region_dist = select_all(regions_coverage.region_dist)
        Array[File] mosdepth_region_summary = regions_coverage.summary
        File multiqc_report = multiqc.multiqc_report
        File? orig_read_count = subsample.orig_read_count
        File? kraken_sequences = kraken.sequences
        File? mosdepth_dups_marked_global_dist = markdups_post.mosdepth_global_dist
        File? mosdepth_dups_marked_global_summary = markdups_post.mosdepth_global_summary
        Array[File]? mosdepth_dups_marked_region_summary
            = markdups_post.mosdepth_region_summary
        Array[File?]? mosdepth_dups_marked_region_dist
            = markdups_post.mosdepth_region_dist
        File? inferred_strandedness = strandedness.strandedness_file
        File? qualimap_rnaseq_results = qualimap_rnaseq.results
        File? junction_summary = junction_annotation.junction_summary
        File? junctions = junction_annotation.junctions
        IntermediateFiles? intermediate_files = optional_files
    }
}

task parse_input {
    meta {
        description: "Parses and validates the `quality_check` workflow's provided inputs"
        outputs: {
            check: "Dummy output to indicate success and to enable call-caching"
            labels: "An array of labels to use on the result coverage files associated with each coverage BED"
        }
    }

    parameter_meta {
        coverage_labels: "An array of equal length to `coverage_beds_len` which determines the prefix label applied to coverage output files. If an empty array is supplied, defaults of `regions1`, `regions2`, etc. will be used."
        input_molecule: "Must be `DNA` or `RNA`"
        gtf_provided: "Was a GTF supplied by the user? Must be `true` if `input_molecule = RNA`."
        coverage_beds_len: "Length of the provided `coverage_beds` array"
    }

    input {
        Array[String] coverage_labels
        String input_molecule
        Boolean gtf_provided
        Int coverage_beds_len
    }

    Int coverage_labels_len = length(coverage_labels)

    command <<<
        EXITCODE=0

        if [ "~{input_molecule}" != "DNA" ] && [ "~{input_molecule}" != "RNA" ]; then
            >&2 echo "molecule input must be 'DNA' or 'RNA'"
            EXITCODE=1
        fi

        if [ "~{input_molecule}" = "RNA" ] && ! ~{gtf_provided}; then
            >&2 echo "Must supply a GTF if molecule = 'RNA'"
            EXITCODE=1
        fi

        if [ "~{coverage_labels_len}" = 0 ]; then
            for (( i=1; i<=~{coverage_beds_len}; i++ )); do
                echo regions$i >> labels.txt
            done
        elif [ "~{coverage_labels_len}" != "~{coverage_beds_len}" ]; then
            >&2 echo "Unequal amount of coverage BEDs and coverage labels."
            >&2 echo "If no labels are provided, generic labels will be created."
            >&2 echo "Otherwise the exact same amount must be supplied."
            EXITCODE=1
        else
            echo "~{sep('\n', coverage_labels)}" >> labels.txt
        fi

        exit $EXITCODE
    >>>

    output {
        String check = "passed"
        Array[String] labels = read_lines("labels.txt")
    }

    runtime {
        memory: "4 GB"
        disk: "10 GB"
        container: 'ghcr.io/stjudecloud/util:1.3.0'
        maxRetries: 1
    }
}

# TODO does this need documentation?
struct IntermediateFiles {
    File? subsampled_bam
    File? subsampled_bam_index
    File? collated_bam
    File? read_one_fastq_gz
    File? read_two_fastq_gz
    File? singleton_reads_fastq_gz
    File? interleaved_reads_fastq_gz  # TODO this will always be empty? Verify and delete
    File? duplicate_marked_bam
    File? duplicate_marked_bam_index
    File? duplicate_marked_bam_md5
}
