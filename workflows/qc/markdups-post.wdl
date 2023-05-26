## # MarkDuplicates Post

version 1.0

import "../../tools/picard.wdl"
import "../../tools/samtools.wdl"
import "../../tools/mosdepth.wdl"

workflow markdups_post {
    meta {
        description: "This workflow runs QC analyses which are impacted by duplicate marking."
        outputs: {
            insert_size_metrics:
                "`*.txt` output file of `picard collectInsertSizeMetrics`"
            insert_size_metrics_pdf:
                "`*.pdf` output file of `picard collectInsertSizeMetrics`"
            flagstat: "`samtools flagstat` report"
            mosdepth_global_summary:
                "Summary of whole genome coverage produced by `mosdepth`"
            mosdepth_global_dist:
                "Distribution of whole genome coverage produced by `mosdepth`"
            mosdepth_region_summary: "Summaries of coverage corresponding to the regions defined by `coverage_beds` input, produced by `mosdepth`"
            mosdepth_region_dist: "Distributions of coverage corresponding to the regions defined by `coverage_beds` input, produced by `mosdepth`"
        }
    }

    parameter_meta {
        markdups_bam: "Input BAM format file to quality check. Duplicates being marked is not necessary for a succesfull run of this workflow."
        markdups_bam_index: "BAM index file corresponding to the input BAM"
        max_retries: "Number of times to retry failed steps. Overrides task level defaults."
        coverage_beds: "An array of 3 column BEDs which are passed to the `-b` flag of mosdepth, in order to restrict coverage analysis to select regions"
        coverage_labels: "An array of equal length to `coverage_beds` which determines the prefix label applied to the output files. If omitted, defaults of `regions1`, `regions2`, etc. will be used."
        prefix: "Prefix for all results files"
    }

    input {
        File markdups_bam
        File markdups_bam_index
        Int? max_retries
        Array[File] coverage_beds = []
        Array[String] coverage_labels = []
        String prefix = basename(markdups_bam, ".bam")
    }

    call picard.collect_insert_size_metrics { input:
            bam=markdups_bam,
            max_retries=max_retries
        }
    call samtools.flagstat as samtools_flagstat { input:
        bam=markdups_bam,
        max_retries=max_retries
    }

    call mosdepth.coverage as wg_coverage {
        input:
            bam=markdups_bam,
            bam_index=markdups_bam_index,
            prefix=prefix + "." + "whole_genome",
            max_retries=max_retries
    }
    scatter(coverage_pair in zip(coverage_beds, coverage_labels)) {
        call mosdepth.coverage as regions_coverage {
            input:
                bam=markdups_bam,
                bam_index=markdups_bam_index,
                prefix=prefix + "." + coverage_pair.right,
                max_retries=max_retries
        }
    }

    output {
        File insert_size_metrics = collect_insert_size_metrics.insert_size_metrics
        File insert_size_metrics_pdf
            = collect_insert_size_metrics.insert_size_metrics_pdf
        File flagstat = samtools_flagstat.flagstat_report
        File mosdepth_global_summary = wg_coverage.summary
        File mosdepth_global_dist = wg_coverage.global_dist
        Array[File] mosdepth_region_summary = regions_coverage.summary
        Array[File?] mosdepth_region_dist = regions_coverage.region_dist
    }
}