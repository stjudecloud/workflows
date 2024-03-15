# SPDX-License-Identifier: MIT
# Copyright St. Jude Children's Research Hospital
version 1.1

import "../../data_structures/flag_filter.wdl"
import "../../tools/fastqc.wdl" as fastqc_tasks
import "../../tools/fq.wdl"
import "../../tools/kraken2.wdl"
import "../../tools/librarian.wdl" as libraran_tasks
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
        help: "Assumes that input BAM is position-sorted."
        external_help: "https://multiqc.info/"
        outputs: {
            bam_checksum: "STDOUT of the `md5sum` command run on the input BAM that has been redirected to a file",
            validate_sam_file: "Validation report produced by `picard ValidateSamFile`. Validation warnings and errors are logged.",
            flagstat_report: "`samtools flagstat` STDOUT redirected to a file. If `mark_duplicates` is `true`, then this result will be generated from the duplicate marked BAM.",
            fastqc_results: "A gzipped tar archive of all FastQC output files",
            instrument_file: "TSV file containing the `ngsderive isntrument` report for the input BAM file",
            read_length_file: "TSV file containing the `ngsderive readlen` report for the input BAM file",
            inferred_encoding: "TSV file containing the `ngsderive encoding` report for the input BAM file",
            alignment_metrics: {
                description: "The text file output of `picard CollectAlignmentSummaryMetrics`",
                external_help: "http://broadinstitute.github.io/picard/picard-metric-definitions.html#AlignmentSummaryMetrics"
            },
            alignment_metrics_pdf: "The PDF file output of `picard CollectAlignmentSummaryMetrics`",
            insert_size_metrics: {
                description: "The text file output of `picard CollectInsertSizeMetrics`. If `mark_duplicates` is `true`, then this result will be generated from the duplicate marked BAM.",
                external_help: "http://broadinstitute.github.io/picard/picard-metric-definitions.html#InsertSizeMetrics"
            },
            insert_size_metrics_pdf: "The PDF file output of `picard CollectInsertSizeMetrics`. If `mark_duplicates` is `true`, then this result will be generated from the duplicate marked BAM.",
            quality_score_distribution_txt: "The text file output of `picard QualityScoreDistribution`",
            quality_score_distribution_pdf: "The PDF file output of `picard QualityScoreDistribution`",
            phred_scores: "Headered TSV file containing PHRED score statistics",
            kraken_report: {
                description: "A Kraken2 summary report",
                external_help: "https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#sample-report-output-format"
            },
            mosdepth_global_dist: "The `$prefix.mosdepth.global.dist.txt` file contains a cumulative distribution indicating the proportion of total bases that were covered for at least a given coverage value. It does this for each chromosome, and for the whole genome.",
            mosdepth_global_summary: "A summary of mean depths per chromosome",
            mosdepth_region_dist: "The `$prefix.mosdepth.region.dist.txt` file contains a cumulative distribution indicating the proportion of total bases in the region(s) defined by the `coverage_bed` that were covered for at least a given coverage value. There will be one file in this array for each `coverage_beds` input file.",
            mosdepth_region_summary: "A summary of mean depths per chromosome and within specified regions per chromosome. There will be one file in this array for each `coverage_beds` input file.",
            multiqc_report: "A gzipped tar archive of all MultiQC output files",
            orig_read_count: "A TSV report containing the original read count before subsampling",
            kraken_sequences: {
                description: "Detailed Kraken2 output that has been gzipped",
                external_help: "https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#standard-kraken-output-format"
            },
            comparative_kraken_report: "Kraken2 summary report for only the alternatively filtered reads",
            comparative_kraken_sequences: "Detailed Kraken2 output for only the alternatively filtered reads",
            mosdepth_dups_marked_global_dist: "The `$prefix.mosdepth.global.dist.txt` file contains a cumulative distribution indicating the proportion of total bases that were covered for at least a given coverage value. It does this for each chromosome, and for the whole genome. This file is produced from analyzing the duplicate marked BAM (only present if `mark_duplicates = true`).",
            mosdepth_dups_marked_global_summary: "A summary of mean depths per chromosome. This file is produced from analyzing the duplicate marked BAM (only present if `mark_duplicates = true`).",
            mosdepth_dups_marked_region_dist: "The `$prefix.mosdepth.region.dist.txt` file contains a cumulative distribution indicating the proportion of total bases in the region(s) defined by the `coverage_bed` that were covered for at least a given coverage value. There will be one file in this array for each `coverage_beds` input file. This file is produced from analyzing the duplicate marked BAM (only present if `mark_duplicates = true`).",
            mosdepth_dups_marked_region_summary: "A summary of mean depths per chromosome and within specified regions per chromosome. There will be one file in this array for each `coverage_beds` input file. This file is produced from analyzing the duplicate marked BAM (only present if `mark_duplicates = true`).",
            inferred_strandedness: "TSV file containing the `ngsderive strandedness` report",
            qualimap_rnaseq_results: "Gzipped tar archive of all QualiMap output files",
            junction_summary: "TSV file containing the `ngsderive junction-annotation` summary",
            junctions: "TSV file containing a detailed list of annotated junctions",
            librarian_report: "A tar archive containing the `librarian` report and raw data.",
            IntermediateFiles: "Any and all files produced as intermediate during pipeline processing. Only output if `output_intermediate_files = true`."
        }
        allowNestedInputs: true
    }

    parameter_meta {
        bam: "Input BAM format file to quality check"
        bam_index: "BAM index file corresponding to the input BAM"
        kraken_db: "Kraken2 database. Can be generated with `../reference/make-qc-reference.wdl`. Must be a tarball without a root directory."
        standard_filter: "Filter to apply to the input BAM while converting to FASTQ, before running Kraken2 and `librarian` (if `run_librarian == true`). This is a `FlagFilter` object (see ../../data_structures/flag_filter.wdl for more information). By default, it will **remove secondary and supplementary reads** from the created FASTQs. **WARNING:** These filters can be tricky to configure; please read documentation thoroughly before changing the defaults. **WARNING:** If you have set `run_librarian` to `true`, we **strongly** recommend leaving this filter at the default value. `librarian` is trained on a specific set of reads, and changing this filter may produce nonsensical results."
        comparative_filter: "Filter to apply to the input BAM while performing a second FASTQ conversion, before running Kraken2 another time. This is a `FlagFilter` object (see ../../data_structures/flag_filter.wdl for more information). By default, it will **remove unmapped, secondary, and supplementary reads** from the created FASTQs. **WARNING** These filters can be tricky to configure; please read documentation thoroughly before changing the defaults."
        gtf: "GTF features file. Gzipped or uncompressed. **Required** for RNA-Seq data."
        multiqc_config: "YAML file for configuring MultiQC"
        extra_multiqc_inputs: "An array of additional files to pass directly into MultiQC"
        coverage_beds: "An array of 3 column BEDs which are passed to the `-b` flag of mosdepth, in order to restrict coverage analysis to select regions. Any regional analysis enabled by this option is _in addition_ to whole genome coverage, which is calculated regardless of this setting. An exon BED and a Coding Sequence BED are examples of regions you may wish to restrict coverage analysis to. Those two BEDs can be created with the workflow in `../reference/make-qc-reference.wdl`."
        coverage_labels: "An array of equal length to `coverage_beds` which determines the prefix label applied to the output files. If omitted, defaults of `regions1`, `regions2`, etc. will be used. If using the BEDs created by `../reference/make-qc-reference.wdl`, the labels [\"exon\", \"CDS\"] are appropriate. Make sure to provide the coverage BEDs **in the same order** as the labels."
        prefix: "Prefix for all results files"
        rna: "Is the sequenced molecule RNA? Enabling this option adds RNA-Seq specific analyses to the workflow. If `true`, a GTF file must be provided. If `false`, the GTF file is ignored."
        mark_duplicates: "Mark duplicates before select analyses? Default behavior is to set this to the value of the `rna` parameter. This is because DNA files are often duplicate marked already, and RNA-Seq files are usually _not_ duplicate marked. If set to `true`, a BAM will be generated and passed to selected downstream analyses. For more details about what analyses are run, review `./markdups-post.wdl`. **WARNING, this duplicate marked BAM is _not_ ouput by default.** If you would like to output this file, set `output_intermediate_files = true`."
        run_librarian: {
            description: "Run the `librarian` tool to generate a report of the likely Illumina library prep kit used to generate the data. **WARNING** this tool is not guaranteed to work on all data, and may produce nonsensical results. `librarian` was trained on a limited set of GEO read data (Gene Expression Oriented). This means the input data should be Paired-End, of mouse or human origin, read length should be >50bp, and derived from a library prep kit that is in the `librarian` database. By default, this tool is run when `rna == true`.",
            external_help: "https://f1000research.com/articles/11-1122/v2",
        }
        run_comparative_kraken: "Run Kraken2 a second time with different FASTQ filtering? If `true`, `comparative_filter` is used in a second run of BAM->FASTQ conversion, resulting in differently filtered FASTQs analyzed by Kraken2. If `false`, `comparative_filter` is ignored."
        output_intermediate_files: "Output intermediate files? FASTQs; if `rna == true` a collated BAM; if `mark_duplicates == true` a duplicate marked BAM, various accessory files like indexes and md5sums; if subsampling was requested _and_ performed then a sampled BAM and associated index. **WARNING, these files can be large.**"
        use_all_cores: "Use all cores? Recommended for cloud environments."
        subsample_n_reads: "Only process a random sampling of approximately `n` reads. Any `n <= 0` for processing entire input. Subsampling is done probabalistically so the exact number of reads in the output will have some variation."
    }

    input {
        File bam
        File bam_index
        File kraken_db
        File? gtf
        File multiqc_config
            = "https://raw.githubusercontent.com/stjudecloud/workflows/main/workflows/qc/inputs/multiqc_config_hg38.yaml"
        Array[File] extra_multiqc_inputs = []
        Array[File] coverage_beds = []
        Array[String] coverage_labels = []
        FlagFilter standard_filter = {
            "include_if_all": "0x0",
            "exclude_if_any": "0x900",  # 0x100 (secondary) || 0x800 (supplementary)
            "include_if_any": "0x0",
            "exclude_if_all": "0x0",
        }
        # TODO: consider making this an array of FlagFilters similar to coverage_beds
        FlagFilter comparative_filter =  {
            "include_if_all": "0x0",
            # 0x4 (unmapped) || 0x100 (secondary) || 0x800 (supplementary)
            "exclude_if_any": "0x904",
            "include_if_any": "0x0",
            "exclude_if_all": "0x0",
        }
        String prefix = basename(bam, ".bam")
        Boolean rna = false
        Boolean mark_duplicates = rna
        Boolean run_librarian = rna
        Boolean run_comparative_kraken = false
        Boolean output_intermediate_files = false
        Boolean use_all_cores = false
        Int subsample_n_reads = -1
    }

    call parse_input { input:
        gtf_provided=defined(gtf),
        rna,
        coverage_beds_len=length(coverage_beds),
        coverage_labels=coverage_labels,
    }
    call flag_filter.validate_FlagFilter as kraken_filter_validator { input:
        flags = standard_filter
    }
    if (run_comparative_kraken) {
        call flag_filter.validate_FlagFilter as comparative_kraken_filter_validator { input:
            flags = comparative_filter
        }
    }

    call md5sum.compute_checksum after parse_input { input: file=bam }

    call samtools.quickcheck after parse_input { input: bam=bam }
    call util.compression_integrity after parse_input { input: bgzipped_file=bam }

    if (subsample_n_reads > 0) {
        call samtools.subsample after quickcheck { input:
            bam=bam,
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
    # If subsampling is disabled **or** input BAM has fewer reads than
    # `subsample_n_reads` this will be `bam`
    File post_subsample_bam = select_first([
        subsample.sampled_bam,
        bam
    ])
    File post_subsample_bam_index = select_first([
        subsample_index.bam_index,
        bam_index
    ])
    String post_subsample_prefix = if (defined(subsample.sampled_bam))
        then prefix + ".subsampled"
        else prefix

    call picard.validate_bam after quickcheck { input:
        bam=post_subsample_bam,
        outfile_name=post_subsample_prefix + ".ValidateSamFile.txt",
        succeed_on_errors=true,
        ignore_list=[],
        summary_mode=true,
    }

    call picard.collect_alignment_summary_metrics after quickcheck { input:
        bam=post_subsample_bam,
        prefix=post_subsample_prefix + ".CollectAlignmentSummaryMetrics",
    }
    call picard.quality_score_distribution after quickcheck { input:
        bam=post_subsample_bam,
        prefix=post_subsample_prefix + ".QualityScoreDistribution",
    }
    call fastqc_tasks.fastqc after quickcheck { input:
        bam=post_subsample_bam,
        prefix=post_subsample_prefix + ".fastqc_results",
        use_all_cores=use_all_cores,
    }
    call ngsderive.instrument after quickcheck { input:
        bam=post_subsample_bam,
        outfile_name=post_subsample_prefix + ".instrument.tsv",
    }
    call ngsderive.read_length after quickcheck { input:
        bam=post_subsample_bam,
        bam_index=post_subsample_bam_index,
        outfile_name=post_subsample_prefix + ".readlength.tsv",
    }
    call ngsderive.encoding after quickcheck { input:
        ngs_files=[post_subsample_bam],
        outfile_name=post_subsample_prefix + ".encoding.tsv",
        num_reads=-1,
    }
    call ngsderive.endedness after quickcheck { input:
        bam=post_subsample_bam,
        outfile_name=post_subsample_prefix + ".endedness.tsv",
        lenient=true,
    }
    call util.global_phred_scores after quickcheck { input:
        bam=post_subsample_bam,
        prefix=post_subsample_prefix,
    }

    call samtools.bam_to_fastq after quickcheck
        after kraken_filter_validator
    { input:
        bam = post_subsample_bam,
        bitwise_filter = standard_filter,
        prefix = post_subsample_prefix,
        # RNA needs a collated BAM for Qualimap
        # DNA can skip the associated storage costs
        retain_collated_bam = rna,
        # disabling fast_mode enables writing of secondary and supplementary alignments
        # to the collated BAM when processing RNA.
        # Those alignments are used downstream by Qualimap.
        fast_mode = (!rna),
        paired_end = true,  # matches default but prevents user from overriding
        interleaved = false,  # matches default but prevents user from overriding
        use_all_cores = use_all_cores,
    }

    call fq.fqlint { input:
        read_one_fastq = select_first([bam_to_fastq.read_one_fastq_gz, "undefined"]),
        read_two_fastq = select_first([bam_to_fastq.read_two_fastq_gz, "undefined"]),
    }
    call kraken2.kraken after fqlint { input:
        read_one_fastq_gz
            = select_first([bam_to_fastq.read_one_fastq_gz, "undefined"]),
        read_two_fastq_gz
            = select_first([bam_to_fastq.read_two_fastq_gz, "undefined"]),
        db=kraken_db,
        prefix=post_subsample_prefix,
        use_all_cores=use_all_cores,
    }
    if (run_librarian) {
        call libraran_tasks.librarian after fqlint { input:
            read_one_fastq = select_first([bam_to_fastq.read_one_fastq_gz, "undefined"]),
        }
    }

    if (run_comparative_kraken) {
        call samtools.bam_to_fastq as alt_filtered_fastq after quickcheck
            after comparative_kraken_filter_validator
        { input:
            bam = post_subsample_bam,
            bitwise_filter = comparative_filter,
            prefix = post_subsample_prefix + ".alt_filtered",
            # matches default but prevents user from overriding
            # If the user wants a collated BAM, they should save the one
            # from the first bam_to_fastq call.
            retain_collated_bam = false,
            # matches default but prevents user from overriding
            # Since the only output here is FASTQs, we can disable fast mode.
            # This discards secondary and supplementary alignments, which should not
            # be converted to FASTQs. (Is that true?)
            fast_mode = true,
            paired_end = true,  # matches default but prevents user from overriding
            interleaved = false,  # matches default but prevents user from overriding
            use_all_cores = use_all_cores,
        }
        call fq.fqlint as alt_filtered_fqlint { input:
            read_one_fastq
                = select_first([alt_filtered_fastq.read_one_fastq_gz, "undefined"]),
            read_two_fastq
                = select_first([alt_filtered_fastq.read_two_fastq_gz, "undefined"]),
        }
        call kraken2.kraken as comparative_kraken after alt_filtered_fqlint { input:
            read_one_fastq_gz
                = select_first([alt_filtered_fastq.read_one_fastq_gz, "undefined"]),
            read_two_fastq_gz
                = select_first([alt_filtered_fastq.read_two_fastq_gz, "undefined"]),
            db = kraken_db,
            prefix = post_subsample_prefix + ".alt_filtered",
            use_all_cores = use_all_cores,
        }
    }

    call mosdepth.coverage as wg_coverage after quickcheck { input:
        bam=post_subsample_bam,
        bam_index=post_subsample_bam_index,
        prefix=post_subsample_prefix + ".whole_genome",
    }
    scatter(coverage_pair in zip(coverage_beds, parse_input.labels)) {
        call mosdepth.coverage as regions_coverage after quickcheck  { input:
            bam=post_subsample_bam,
            bam_index=post_subsample_bam_index,
            coverage_bed=coverage_pair.left,
            prefix=post_subsample_prefix + "." + coverage_pair.right,
        }
    }

    if (rna) {
        call ngsderive.junction_annotation after quickcheck { input:
            bam=post_subsample_bam,
            bam_index=post_subsample_bam_index,
            gene_model=select_first([gtf, "undefined"]),
            prefix=post_subsample_prefix,
        }
        call ngsderive.strandedness after quickcheck { input:
            bam=post_subsample_bam,
            bam_index=post_subsample_bam_index,
            gene_model=select_first([gtf, "undefined"]),
            outfile_name=post_subsample_prefix + ".strandedness.tsv",
        }
        call qualimap.rnaseq as qualimap_rnaseq { input:
            bam=select_first([bam_to_fastq.collated_bam, "undefined"]),
            prefix=post_subsample_prefix + ".qualimap_rnaseq_results",
            gtf=select_first([gtf, "undefined"]),
            name_sorted=true,
            paired_end=true,  # matches default but prevents user from overriding
        }
    }
    if (mark_duplicates) {
        call picard.mark_duplicates as markdups after quickcheck { input:
            bam = post_subsample_bam,
            create_bam = true,
            prefix = post_subsample_prefix + ".MarkDuplicates",
            optical_distance = 0,
        }
        call markdups_post_wf.markdups_post { input:
            markdups_bam = select_first([
                markdups.duplicate_marked_bam,
                "undefined",
            ]),
            markdups_bam_index = select_first([
                markdups.duplicate_marked_bam_index,
                "undefined",
            ]),
            coverage_beds = coverage_beds,
            coverage_labels = parse_input.labels,
            prefix = post_subsample_prefix + ".MarkDuplicates",
        }
    }
    if (!mark_duplicates) {
        # These analyses are called in the markdups_post workflow.
        # They should still be run if duplicates were not marked.
        call picard.collect_insert_size_metrics after quickcheck { input:
            bam = post_subsample_bam,
            prefix = post_subsample_prefix + ".CollectInsertSizeMetrics",
        }
        call samtools.flagstat after quickcheck { input:
            bam = post_subsample_bam,
            outfile_name = post_subsample_prefix + ".flagstat.txt",
        }
    }

    call multiqc_tasks.multiqc { input:
        input_files=select_all(flatten([
            [
                validate_bam.validate_report,
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
                librarian.raw_data,
                comparative_kraken.report,
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
            "sampled_bam": subsample.sampled_bam,
            "sampled_bam_index": subsample_index.bam_index,
            "collated_bam": bam_to_fastq.collated_bam,
            "read_one_fastq_gz": bam_to_fastq.read_one_fastq_gz,
            "read_two_fastq_gz": bam_to_fastq.read_two_fastq_gz,
            "singleton_reads_fastq_gz": bam_to_fastq.singleton_reads_fastq_gz,
            "alt_filtered_read_one_fastq_gz": alt_filtered_fastq.read_one_fastq_gz,
            "alt_filtered_read_two_fastq_gz": alt_filtered_fastq.read_two_fastq_gz,
            "alt_filtered_singleton_reads_fastq_gz": alt_filtered_fastq.singleton_reads_fastq_gz,
            "duplicate_marked_bam": markdups.duplicate_marked_bam,
            "duplicate_marked_bam_index": markdups.duplicate_marked_bam_index,
            "duplicate_marked_bam_md5": markdups.duplicate_marked_bam_md5,
        }
    }

    output {
        File bam_checksum = compute_checksum.md5sum
        File validate_sam_file = validate_bam.validate_report
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
        File? comparative_kraken_report = comparative_kraken.report
        File? comparative_kraken_sequences = comparative_kraken.sequences
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
        File? librarian_report = librarian.report
        IntermediateFiles? intermediate_files = optional_files
    }
}

task parse_input {
    meta {
        description: "Parses and validates the `quality_check` workflow's provided inputs"
        outputs: {
            check: "Dummy output to indicate success and to enable call-caching",  # TODO: is this still needed with labels?
            labels: "An array of labels to use on the result coverage files associated with each coverage BED"
        }
    }

    parameter_meta {
        coverage_labels: "An array of equal length to `coverage_beds_len` which determines the prefix label applied to coverage output files. If an empty array is supplied, defaults of `regions1`, `regions2`, etc. will be used."
        rna: "Is the sequenced molecule RNA?"
        gtf_provided: "Was a GTF supplied by the user? Must be `true` if `rna == true`."
        coverage_beds_len: "Length of the provided `coverage_beds` array"
    }

    input {
        Array[String] coverage_labels
        Boolean rna
        Boolean gtf_provided
        Int coverage_beds_len
    }

    Int coverage_labels_len = length(coverage_labels)

    command <<<
        EXITCODE=0

        if ~{rna} && ! ~{gtf_provided}; then
            >&2 echo "Must supply a GTF if 'rna == true'"
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
    File? sampled_bam
    File? sampled_bam_index
    File? collated_bam
    File? read_one_fastq_gz
    File? read_two_fastq_gz
    File? singleton_reads_fastq_gz
    File? alt_filtered_read_one_fastq_gz
    File? alt_filtered_read_two_fastq_gz
    File? alt_filtered_singleton_reads_fastq_gz
    File? duplicate_marked_bam
    File? duplicate_marked_bam_index
    File? duplicate_marked_bam_md5
}
