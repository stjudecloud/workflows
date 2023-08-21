## # Quality Check Standard
##
## This workflow runs a variety of quality checking software on any BAM file.
## It can be WGS, WES, or Transcriptome data. The results are aggregated and
## run through [MultiQC](https://multiqc.info/).
##
## ## LICENSING
## 
## #### MIT License
##
## Copyright 2020-Present St. Jude Children's Research Hospital
##
## Permission is hereby granted, free of charge, to any person obtaining a copy of this
## software and associated documentation files (the "Software"), to deal in the Software
## without restriction, including without limitation the rights to use, copy, modify, merge,
## publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons
## to whom the Software is furnished to do so, subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in all copies or
## substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
## BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
## NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
## DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
## OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

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
    parameter_meta {
        bam: "Input BAM format file to quality check"
        bam_index: "BAM index file corresponding to the input BAM"
        reference_fasta: "Reference genome in FASTA format"
        kraken_db: "Kraken2 database. Can be generated with `make-qc-reference.wdl`. Must be a tarball without a root directory."
        molecule: {
            description: "Data type"
            choices: [
                'DNA',
                'RNA'
            ]
        }
        gtf: "GTF features file. Gzipped or uncompressed. **Required** for RNA-Seq data."
        max_retries: "Number of times to retry failed steps. Overrides task level defaults."
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
        File reference_fasta
        File kraken_db
        String molecule
        File? gtf
        Int? max_retries
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
        max_retries=max_retries
    }

    call md5sum.compute_checksum { input: file=bam, max_retries=max_retries }

    call samtools.quickcheck { input: bam=bam, max_retries=max_retries }
    call util.compression_integrity { input: bam=bam, max_retries=max_retries }

    if (subsample_n_reads > 0) {
        call samtools.subsample { input:
            bam=quickcheck.checked_bam,
            prefix=prefix,
            desired_reads=subsample_n_reads,
            use_all_cores=use_all_cores,
            max_retries=max_retries
        }
        if (defined(subsample.sampled_bam)) {
            call samtools.index as subsample_index { input:
                bam=select_first([subsample.sampled_bam, "undefined"]),
                use_all_cores=use_all_cores,
                max_retries=max_retries
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
        reference_fasta=reference_fasta,
        succeed_on_errors=true,
        ignore_list=[],
        summary_mode=true,
        max_retries=max_retries
    }

    call picard.collect_alignment_summary_metrics { input:
        bam=post_subsample_bam,
        prefix=post_subsample_prefix + ".CollectAlignmentSummaryMetrics",
        max_retries=max_retries
    }
    call picard.collect_gc_bias_metrics { input:
        bam=post_subsample_bam,
        reference_fasta=reference_fasta,
        prefix=post_subsample_prefix + ".CollectGcBiasMetrics",
        max_retries=max_retries
    }
    call picard.quality_score_distribution { input:
        bam=post_subsample_bam,
        prefix=post_subsample_prefix + ".QualityScoreDistribution",
        max_retries=max_retries
    }
    call fastqc_tasks.fastqc { input:
        bam=post_subsample_bam,
        prefix=post_subsample_prefix + ".fastqc_results",
        use_all_cores=use_all_cores,
        max_retries=max_retries
    }
    call ngsderive.instrument { input:
        bam=post_subsample_bam,
        outfile_name=post_subsample_prefix + ".instrument.tsv",
        max_retries=max_retries
    }
    call ngsderive.read_length { input:
        bam=post_subsample_bam,
        bam_index=post_subsample_bam_index,
        outfile_name=post_subsample_prefix + ".readlength.tsv",
        max_retries=max_retries
    }
    call ngsderive.encoding { input:
        ngs_files=[post_subsample_bam],
        outfile_name=post_subsample_prefix + ".encoding.tsv",
        max_retries=max_retries
    }
    call util.global_phred_scores { input:
        bam=post_subsample_bam,
        prefix=post_subsample_prefix,
        max_retries=max_retries
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
        max_retries=max_retries
    }

    call fq.fqlint { input:
        read_one_fastq=select_first([collate_to_fastq.read_one_fastq_gz, "undefined"]),
        read_two_fastq=collate_to_fastq.read_two_fastq_gz,
        max_retries=max_retries
    }
    call kraken2.kraken { input:
        read_one_fastq_gz=fqlint.validated_read1,
        read_two_fastq_gz=select_first([fqlint.validated_read2, "undefined"]),
        db=kraken_db,
        prefix=post_subsample_prefix,
        use_all_cores=use_all_cores,
        max_retries=max_retries
    }

    call mosdepth.coverage as wg_coverage { input:
        bam=post_subsample_bam,
        bam_index=post_subsample_bam_index,
        prefix=post_subsample_prefix + ".whole_genome",
        max_retries=max_retries
    }
    scatter(coverage_pair in zip(coverage_beds, parse_input.labels)) {
        call mosdepth.coverage as regions_coverage { input:
            bam=post_subsample_bam,
            bam_index=post_subsample_bam_index,
            coverage_bed=coverage_pair.left,
            prefix=post_subsample_prefix + "." + coverage_pair.right,
            max_retries=max_retries
        }
    }

    if (molecule == "RNA") {
        call ngsderive.junction_annotation { input:
            bam=post_subsample_bam,
            bam_index=post_subsample_bam_index,
            gtf=select_first([gtf, "undefined"]),
            prefix=post_subsample_prefix,
            max_retries=max_retries
        }
        call ngsderive.infer_strandedness { input:
            bam=post_subsample_bam,
            bam_index=post_subsample_bam_index,
            gtf=select_first([gtf, "undefined"]),
            outfile_name=post_subsample_prefix + ".strandedness.tsv",
            max_retries=max_retries
        }
        call qualimap.rnaseq as qualimap_rnaseq { input:
            bam=select_first([collate_to_fastq.collated_bam, "undefined"]),
            prefix=post_subsample_prefix + ".qualimap_rnaseq_results",
            gtf=select_first([gtf, "undefined"]),
            name_sorted=true,
            paired_end=true,  # matches default but prevents user from overriding
            max_retries=max_retries
        }
    }

    call picard.mark_duplicates as markdups { input:
        bam=post_subsample_bam,
        create_bam=mark_duplicates,
        prefix=post_subsample_prefix + ".MarkDuplicates",
        max_retries=max_retries
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
            max_retries=max_retries
        }
    }
    if (! mark_duplicates) {
        # These analyses are called in the markdups_post workflow.
        # They should still be run if duplicates were not marked.
        call picard.collect_insert_size_metrics { input:
            bam=post_subsample_bam,
            prefix=post_subsample_prefix + ".CollectInsertSizeMetrics",
            max_retries=max_retries
        }
        call samtools.flagstat { input:
            bam=post_subsample_bam,
            outfile_name=post_subsample_prefix + ".flagstat.txt",
            max_retries=max_retries
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
                fastqc.raw_data,
                collect_alignment_summary_metrics.alignment_metrics,
                collect_gc_bias_metrics.gc_bias_metrics,
                collect_gc_bias_metrics.gc_bias_metrics_summary,
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
                infer_strandedness.strandedness_file,
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
        max_retries=max_retries
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
        File alignment_metrics = collect_alignment_summary_metrics.alignment_metrics
        File alignment_metrics_pdf
            = collect_alignment_summary_metrics.alignment_metrics_pdf
        File gc_bias_metrics = collect_gc_bias_metrics.gc_bias_metrics
        File gc_bias_metrics_summary = collect_gc_bias_metrics.gc_bias_metrics_summary
        File gc_bias_metrics_pdf = collect_gc_bias_metrics.gc_bias_metrics_pdf
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
        File? inferred_strandedness = infer_strandedness.strandedness_file
        File? qualimap_rnaseq_results = qualimap_rnaseq.results
        File? junction_summary = junction_annotation.junction_summary
        File? junctions = junction_annotation.junctions
        IntermediateFiles? intermediate_files = optional_files
    }
}

task parse_input {
    input {
        Array[String] coverage_labels
        String input_molecule
        Boolean gtf_provided
        Int coverage_beds_len
        Int memory_gb = 4
        Int disk_size_gb = 10
        Int max_retries = 1
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
        memory: "~{memory_gb} GB"
        disk: "~{disk_size_gb} GB"
        docker: 'ghcr.io/stjudecloud/util:1.3.0'
        maxRetries: max_retries
    }
}

struct IntermediateFiles {
    File? subsampled_bam
    File? subsampled_bam_index
    File? collated_bam
    File? read_one_fastq_gz
    File? read_two_fastq_gz
    File? singleton_reads_fastq_gz
    File? interleaved_reads_fastq_gz
    File? duplicate_marked_bam
    File? duplicate_marked_bam_index
    File? duplicate_marked_bam_md5
}
