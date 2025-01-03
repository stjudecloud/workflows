version 1.1

import "../../tools/hilow.wdl"
import "../../tools/samtools.wdl"

workflow hic_post {
    meta {
        description: "Post-processing of Hi-C data."
        outputs: {
            hic_file: "Interactions in .hic format",
            filtered_pairs: "File with filtered list of interaction pairs",
            qc_report: "QC report for Hi-C experiment",
            combined_bam: "Merged and duplicate marked BAM file",
            combined_bam_index: "Index for the combined BAM file",
        }
    }

    parameter_meta {
        all_valid_pairs: "All valid pairs file from HiC-Pro"
        chromsizes: {
            description: "Tab delimited file with chromosome sizes"
        }
        exclude_list: "BED file with regions to exclude from analysis"
        prefix: "Prefix for the output BAM"
        contact_stats: "Contact stats file"
        merged_read_one_mapping_stats: "Mapping stats for read one from HiC-Pro"
        merged_read_two_mapping_stats: "Mapping stats for read two from HiC-Pro"
        merged_pairing_stats: "Pairing stats for merged reads from HiC-Pro"
        combined_bams: "Array of BAM files to merge"
    }

    input {
        File all_valid_pairs
        File chromsizes
        File contact_stats
        File merged_read_one_mapping_stats
        File merged_read_two_mapping_stats
        File merged_pairing_stats
        Array[File] combined_bams
        String prefix
        File? exclude_list
    }

    call hilow.converthic { input:
        all_valid_pairs = all_valid_pairs,
        chromsizes,
    }

    if (defined(exclude_list)) {
        call hilow.filter { input:
            all_valid_pairs = all_valid_pairs,
            chromsizes,
            exclude_list = select_first([exclude_list, ""]),
        }
    }

    call hilow.qcreport { input:
        prefix,
        all_valid_pairs_stats = contact_stats,
        mapping_stats_read1 = merged_read_one_mapping_stats,
        mapping_stats_read2 = merged_read_two_mapping_stats,
        pairing_stats = merged_pairing_stats,
    }

    call samtools.merge { input:
        bams = combined_bams,
        prefix = prefix + ".merged",
        attach_rg = false,
    }

    call samtools.fixmate { input:
        bam = merge.merged_bam,
    }

    call samtools.sort { input:
        bam = fixmate.fixmate_bam,
    }

    call samtools.markdup { input:
        bam = sort.sorted_bam,
        create_bam = true,
        mark_supp_or_sec_or_unmapped_as_duplicates = true,
        prefix,
    }

    call samtools.index { input:
        bam = select_first([markdup.markdup_bam, ""]),
    }

    output {
        File hic_file = converthic.hic_file
        File? filtered_pairs = filter.filtered_pairs
        File qc_report = qcreport.qc_report
        File combined_bam = select_first([markdup.markdup_bam, ""])
        File combined_bam_index = index.bam_index
    }
}
