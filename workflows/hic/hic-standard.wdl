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
            qc_report: "QC report",
            mapping_stats_plot: "Mapping stats plot",
            pairing_stats_plot: "Pairing stats plot",
            filtering_stats_plot: "Filtering stats plot",
            filtering_size_plot: "Filtering size plot",
            contacts_stats_plot: "Contacts stats plot",
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

    # scatter (pair in zip(flatten(r1_split.fastqs), flatten(r2_split.fastqs))) {
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

    call hilow.qc_hic { input:
        plot_type = "all",
        mapping_stats = flatten([
            hicpro_core.read_one_mapping_stats,
            hicpro_core.read_two_mapping_stats,
        ]),
        pairing_stats = hicpro_core.pairing_stats,
        fragment_stats = flatten([
            hicpro_core.rs_stats,
            hicpro_core.fragment_stats,
        ]),
        contacts_stats = [hicpro_core.contact_stats],
        sample_name = prefix,
        remove_singleton = true,
        remove_multimapper = true,
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
        mapping_stats_R1 = hicpro_core.merged_read_one_mapping_stats,
        mapping_stats_R2 = hicpro_core.merged_read_two_mapping_stats,
        pairing_stats = hicpro_core.merged_pairing_stats,
    }

    call samtools.merge { input:
        bams = hicpro_core.combined_bams,
    }

    call samtools.markdup { input:
        bam = merge.merged_bam,
        mark_supp_or_sec_or_unmapped_as_duplicates = true,
    }

    output {
        File hic_file = converthic.hic_file
        File? filtered_pairs = filter.filtered_pairs
        File all_valid_pairs = hicpro_core.all_valid_pairs
        File qc_report = qcreport.qc_report
        File? mapping_stats_plot = qc_hic.mapping_stats_plot
        File? pairing_stats_plot = qc_hic.pairing_stats_plot
        File? filtering_stats_plot = qc_hic.filtering_stats_plot
        File? filtering_size_plot = qc_hic.filtering_size_plot
        File? contacts_stats_plot = qc_hic.contacts_stats_plot
        Array[File] ice_normalized_matrices = hicpro_core.ice_normalized_matrices
        File combined_bam = select_first([markdup.markdup_bam, ""])
    }
}
