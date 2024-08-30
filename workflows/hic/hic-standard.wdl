version 1.1

import "../../data_structures/flag_filter.wdl"
import "../../data_structures/read_group.wdl"
import "../../tools/bowtie2.wdl"
import "../../tools/hilow.wdl"
import "../../tools/samtools.wdl"
import "../../tools/util.wdl"
import "hicpro-core.wdl" as hicpro_core

workflow hic_standard {
    meta {
        description: "Standard Hi-C processing workflow."
        outputs: {
            hic_file: "Hi-C file in HIC format",
            filtered_pairs: "Filtered pairs file",
            all_valid_pairs: "All valid pairs file",
            qc_report: "QC report for Hi-C experiment. Contains total pairs, R1 and R2 alignment, valid pairs, and interactions along with categorical classificaiton.",
            mapping_stats_plot: "Mapping stats plot for R1 and R2 tags. Contains aligned percentage, unaligned percentage, and full-read and trimmed read mapping.",
            pairing_stats_plot: "Pairing stats plot containing reported status for all pairs and quality for filtered pairs.",
            filtering_stats_plot: "Plot of pair alignments to restriction fragments. Contains orientation and classification of intercations.",
            filtering_size_plot: "Plot of distribution of fragment sizes.",
            contacts_stats_plot: "Plot of contact ranges for valid pairs.",
            ice_normalized_matrices: "ICE normalized matrices",
            combined_bam: "Samtools merged and duplicate marked BAM file",
        }
    }

    parameter_meta {
        read_one_fastqs_gz: "An array of gzipped FASTQ files containing read one information"
        bowtie_db_tar_gz: "A gzipped TAR file containing the bowtie2 reference files."
        prefix: "Prefix for the BAM"
        read_two_fastqs_gz: {
            description: "An array of gzipped FASTQ files containing read two information",
            common: true,
        }
        capture_bed: "BED file of target regions for capture Hi-C data"
        allele_specific_snp: {
            description: "VCF file of SNPs to use in distinguishing parental origin",
            external_help: "https://nservant.github.io/HiC-Pro/AS.html#as",
        }
        exclude_list: "BED file with regions to exclude from analysis"
        fragment_file: "BED file with restriction fragments"
        ligation_site: "Ligation site sequence used for reads trimming."
        bin_sizes: {
            description: "Resolution of contact maps to generate ",
            hicpro_field: "BIN_SIZE",
        }
        chromsizes: {
            description: "Tab delimited file with chromosome sizes"
        }
    }

    input {
        File bowtie_db_tar_gz
        File chromsizes
        Array[File] read_one_fastqs_gz
        String prefix
        File? exclude_list
        File? fragment_file
        File? capture_bed
        File? allele_specific_snp
        String? ligation_site = "GATCGATC"
        Array[File] read_two_fastqs_gz = []
        Array[Int] bin_sizes = [
            5000,
            10000,
            25000,
            50000,
            100000,
            250000,
            500000,
            1000000,
            2500000,
        ]
    }

    scatter (pair in zip(read_one_fastqs_gz, read_two_fastqs_gz)) {
        call util.split_fastq as r1_split { input:
            fastq = pair.left,
            reads_per_file = 30000000,
        }
        call util.split_fastq as r2_split { input:
            fastq = pair.right,
            reads_per_file = 30000000,
        }
    }

    call hicpro_core.hicpro_core { input:
        read_one_fastqs_gz = flatten(r1_split.fastqs),
        read_two_fastqs_gz = flatten(r2_split.fastqs),
        bowtie_db_tar_gz,
        chromsizes,
        fragment_file,
        capture_bed,
        allele_specific_snp,
        prefix,
        ligation_site,
        bin_sizes,
        matrix_format = "upper",
        remove_singleton = true,
        remove_multimapper = true,
        remove_duplicates = true,
    }

    call hilow.converthic { input:
        all_valid_pairs = hicpro_core.all_valid_pairs,
        chromsizes,
    }

    if (defined(exclude_list)) {
        call hilow.filter { input:
            all_valid_pairs = hicpro_core.all_valid_pairs,
            chromsizes,
            exclude_list = select_first([exclude_list, ""]),
        }
    }

    call hilow.qcreport { input:
        prefix,
        all_valid_pairs_stats = hicpro_core.contact_stats,
        mapping_stats_read1 = hicpro_core.merged_read_one_mapping_stats,
        mapping_stats_read2 = hicpro_core.merged_read_two_mapping_stats,
        pairing_stats = hicpro_core.merged_pairing_stats,
    }

    call samtools.merge { input:
        bams = hicpro_core.combined_bams,
        prefix = prefix + ".merged",
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
        File all_valid_pairs = hicpro_core.all_valid_pairs
        File qc_report = qcreport.qc_report
        File? mapping_stats_plot = hicpro_core.mapping_stats_plot
        File? pairing_stats_plot = hicpro_core.pairing_stats_plot
        File? filtering_stats_plot = hicpro_core.filtering_stats_plot
        File? filtering_size_plot = hicpro_core.filtering_size_plot
        File? contacts_stats_plot = hicpro_core.contacts_stats_plot
        Array[File] ice_normalized_matrices = hicpro_core.ice_normalized_matrices
        File combined_bam = select_first([markdup.markdup_bam, ""])
        File combined_bam_index = index.bam_index
    }
}
