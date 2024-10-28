version 1.1

import "../../data_structures/flag_filter.wdl"
import "../../data_structures/read_group.wdl"
import "../../tools/bowtie2.wdl"
import "../../tools/hilow.wdl"
import "../../tools/samtools.wdl"
import "../../tools/util.wdl"
import "hicpro-core.wdl" as hicpro_core
import "hic-post.wdl" as hicpro_post

workflow hic_standard_fastq {
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
            combined_bam_index: "Index for the combined BAM file",
        }
    }

    parameter_meta {
        bowtie_db_tar_gz: "A gzipped TAR file containing the bowtie2 reference files."
        chromsizes: {
            description: "Tab delimited file with chromosome sizes. File is a headerless, two-column TSV and column 1 is the chromosome label and column 2 is the chromosome size in base pairs (bp).",
        }
        read_one_fastqs_gz: "An array of gzipped FASTQ files containing read one information"
        read_groups: "An array of ReadGroup structs containing read group information for each input FASTQ to output in the BAM file"
        prefix: "Prefix for the BAM"
        exclude_list: "BED file with regions to exclude from analysis"
        fragment_file: "BED file with restriction fragments"
        capture_bed: "BED file of target regions for capture Hi-C data"
        allele_specific_snp: {
            description: "VCF file of SNPs to use in distinguishing parental origin",
            external_help: "https://nservant.github.io/HiC-Pro/AS.html#as",
        }
        ligation_site: "Ligation site sequence used for reads trimming."
        read_two_fastqs_gz: {
            description: "An array of gzipped FASTQ files containing read two information",
            common: true,
        }
        bin_sizes: {
            description: "Resolution of contact maps to generate ",
            hicpro_field: "BIN_SIZE",
        }
    }

    input {
        File bowtie_db_tar_gz
        File chromsizes
        Array[File] read_one_fastqs_gz
        Array[ReadGroup] read_groups
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

    call hicpro_core.hicpro_core { input:
        read_one_fastqs_gz = read_one_fastqs_gz,
        read_two_fastqs_gz = read_two_fastqs_gz,
        read_groups,
        bowtie_db_tar_gz,
        chromsizes,
        fragment_file,
        capture_bed,
        allele_specific_snp,
        prefix,
        ligation_site,
        bin_sizes,
    }

    call hicpro_post.hic_post { input:
        all_valid_pairs = hicpro_core.all_valid_pairs,
        chromsizes,
        exclude_list,
        prefix,
        contact_stats = hicpro_core.contact_stats,
        merged_read_one_mapping_stats = hicpro_core.merged_read_one_mapping_stats,
        merged_read_two_mapping_stats = hicpro_core.merged_read_two_mapping_stats,
        merged_pairing_stats = hicpro_core.merged_pairing_stats,
        combined_bams = hicpro_core.combined_bams,
    }

    output {
        File hic_file = hic_post.hic_file
        File? filtered_pairs = hic_post.filtered_pairs
        File all_valid_pairs = hicpro_core.all_valid_pairs
        File qc_report = hic_post.qc_report
        File? mapping_stats_plot = hicpro_core.mapping_stats_plot
        File? pairing_stats_plot = hicpro_core.pairing_stats_plot
        File? filtering_stats_plot = hicpro_core.filtering_stats_plot
        File? filtering_size_plot = hicpro_core.filtering_size_plot
        File? contacts_stats_plot = hicpro_core.contacts_stats_plot
        Array[File] ice_normalized_matrices = hicpro_core.ice_normalized_matrices
        File combined_bam = hic_post.combined_bam
        File combined_bam_index = hic_post.combined_bam_index
    }
}
