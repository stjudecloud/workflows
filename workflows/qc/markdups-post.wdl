## # MarkDuplicates Post
##
## An investigation of all our QC tools was conducted when duplicate marking was
## introduced to our pipeline. Most tools do not take into consideration whether a read
## is a duplicate or not. But the tasks called below produce different results depending
## on whether the input BAM has been duplicate marked or not.
version 1.1

import "../../tools/mosdepth.wdl"
import "../../tools/picard.wdl"
import "../../tools/samtools.wdl"

workflow markdups_post {
    meta {
        name: "Mark Duplicates Post"
        description: "Runs QC analyses which are impacted by duplicate marking"
        outputs: {
            insert_size_metrics: "`*.txt` output file of `picard collectInsertSizeMetrics`",
            insert_size_metrics_pdf: "`*.pdf` output file of `picard collectInsertSizeMetrics`",
            flagstat_report: "`samtools flagstat` report",
            mosdepth_global_summary: "Summary of whole genome coverage produced by `mosdepth`",
            mosdepth_global_dist: "Distribution of whole genome coverage produced by `mosdepth`",
            mosdepth_region_summary: "Summaries of coverage corresponding to the regions defined by `coverage_beds` input, produced by `mosdepth`",
            mosdepth_region_dist: "Distributions of coverage corresponding to the regions defined by `coverage_beds` input, produced by `mosdepth`",
        }
        allowNestedInputs: true
    }

    parameter_meta {
        markdups_bam: "Input BAM format file to quality check. Duplicates being marked is not necessary for a successful run of this workflow."
        markdups_bam_index: "BAM index file corresponding to the input BAM"
        coverage_beds: "An array of 3 column BEDs which are passed to the `-b` flag of mosdepth, in order to restrict coverage analysis to select regions"
        coverage_labels: {
            description: "An array of equal length to `coverage_beds` which determines the prefix label applied to the output files.",
            help: "If omitted, defaults of `regions1`, `regions2`, etc. will be used.",
        }
        prefix: "Prefix for all results files"
    }

    input {
        File markdups_bam
        File markdups_bam_index
        Array[File] coverage_beds = []
        Array[String] coverage_labels = []
        String prefix = basename(markdups_bam, ".bam")
    }

    call picard.collect_insert_size_metrics { input:
        bam = markdups_bam,
        prefix = prefix + ".CollectInsertSizeMetrics",
    }
    call samtools.flagstat { input:
        bam = markdups_bam,
        outfile_name = prefix + ".flagstat.txt",
    }

    call mosdepth.coverage as wg_coverage { input:
        bam = markdups_bam,
        bam_index = markdups_bam_index,
        prefix = prefix + "." + "whole_genome",
    }
    scatter (coverage_pair in zip(coverage_beds, coverage_labels)) {
        call mosdepth.coverage as regions_coverage { input:
            bam = markdups_bam,
            bam_index = markdups_bam_index,
            coverage_bed = coverage_pair.left,
            prefix = prefix + "." + coverage_pair.right,
        }
    }

    output {
        File insert_size_metrics = collect_insert_size_metrics.insert_size_metrics
        File insert_size_metrics_pdf = collect_insert_size_metrics.insert_size_metrics_pdf
        File flagstat_report = flagstat.flagstat_report
        File mosdepth_global_summary = wg_coverage.summary
        File mosdepth_global_dist = wg_coverage.global_dist
        Array[File] mosdepth_region_summary = regions_coverage.summary
        Array[File?] mosdepth_region_dist = regions_coverage.region_dist
    }
}
