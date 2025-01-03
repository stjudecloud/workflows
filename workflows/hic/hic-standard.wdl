version 1.1

import "../../data_structures/read_group.wdl"
import "../../tools/picard.wdl"
import "../../tools/samtools.wdl"
import "../general/bam-to-fastqs.wdl" as bam_to_fastqs_wf
import "hicpro-core.wdl" as hicpro_core
import "hic-post.wdl" as hicpro_post

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
            combined_bam_index: "Index for the combined BAM file",
        }
    }

    parameter_meta {
        bam: "BAM file to extract reads and harmonize"
        bowtie_db_tar_gz: "A gzipped TAR file containing the bowtie2 reference files."
        prefix: "Prefix for the output BAM"
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
        subsample_n_reads: "Only process a random sampling of `n` reads. Any `n`<=`0` for processing entire input."
        validate_input: "Ensure input BAM is well-formed before beginning harmonization?"
    }

    input {
        File bowtie_db_tar_gz
        File chromsizes
        File bam
        String prefix
        File? exclude_list
        File? fragment_file
        File? capture_bed
        File? allele_specific_snp
        String? ligation_site = "GATCGATC"
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
        Boolean validate_input = true
        Int subsample_n_reads = -1
    }

    if (validate_input) {
        call picard.validate_bam as validate_input_bam { input:
            bam = bam,
        }
    }

    if (subsample_n_reads > 0) {
        call samtools.subsample { input:
            bam,
            desired_reads = subsample_n_reads,
        }
    }
    File selected_bam = select_first([subsample.sampled_bam, bam])

    call bam_to_fastqs_wf.bam_to_fastqs { input:
        bam = selected_bam,
    }

    call read_group.get_read_groups { input:
        bam = selected_bam,
    }

    call hicpro_core.hicpro_core { input:
        read_one_fastqs_gz = bam_to_fastqs.read1s,
        read_two_fastqs_gz = select_all(bam_to_fastqs.read2s),
        read_groups = get_read_groups.read_groups,
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
