## # MarkDuplicates Post
##
## This workflow runs QC analyses which are impacted by duplicate marking.

version 1.0

import "../../tools/picard.wdl"
import "../../tools/samtools.wdl"
import "../../tools/mosdepth.wdl"

workflow markdups_post {
    input {
        File markdups_bam
        File markdups_bam_index
        String prefix = basename(markdups_bam, ".bam")
        Array[File] coverage_beds = []
        Array[String] coverage_labels = []
        Int? max_retries
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