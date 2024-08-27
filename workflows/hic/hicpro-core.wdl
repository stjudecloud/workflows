version 1.1

import "../../data_structures/flag_filter.wdl"
import "../../data_structures/read_group.wdl"
import "../../tools/bowtie2.wdl"
import "../../tools/hicpro.wdl"
import "../../tools/hilow.wdl"
import "../../tools/samtools.wdl"
import "../../tools/util.wdl"

workflow hicpro_core {
    meta {
        description: "HiC-Pro implementation."
        outputs: {
            hic_file: "HiC-Pro output file",
            filtered_pairs: "Filtered pairs file",
            all_valid_pairs: "All valid pairs file",
            qc_report: "QC report",
            mapping_stats_plot: "Mapping stats plot",
            pairing_stats_plot: "Pairing stats plot",
            filtering_stats_plot: "Filtering stats plot",
            filtering_size_plot: "Filtering size plot",
            contacts_stats_plot: "Contacts stats plot",
            ice_normalized_matrices: "ICE normalized matrices",
            combined_bams: "Combined BAM files",
        }
    }

    parameter_meta {
        read_one_fastqs_gz: "An array of gzipped FASTQ files containing read one information"
        bowtie_db_tar_gz: "A gzipped TAR file containing the bowtie2 reference files."
        prefix: "Prefix for the BAM"
        read_groups: "A string containing the read group information to output in the BAM file. If including multiple read group fields per-read group, they should be space delimited. Read groups should be comma separated, with a space on each side (i.e. ' , '). The ID field must come first for each read group and must be contained in the basename of a FASTQ file or pair of FASTQ files if Paired-End. Example: `ID:rg1 PU:flowcell1.lane1 SM:sample1 PL:illumina LB:sample1_lib1 , ID:rg2 PU:flowcell1.lane2 SM:sample1 PL:illumina LB:sample1_lib1`. These two read groups could be associated with the following four FASTQs: `sample1.rg1_R1.fastq,sample1.rg2_R1.fastq` and `sample1.rg1_R2.fastq,sample1.rg2_R2.fastq`"
        read_two_fastqs_gz: {
            description: "An array of gzipped FASTQ files containing read two information",
            common: true,
        }
        capture_bed: "BED file of target regions for capture Hi-C data"
        allele_specific_snp: {
            description: "VCF file of SNPs to use in distinguishing parental origin",
            external_help: "https://nservant.github.io/HiC-Pro/AS.html#as",
        }
        genome_prefix: "Reference genome prefix used for genome indexes"
        exclude_list: "BED file with regions to exclude from analysis"
        fragment_file: "BED file with restriction fragments"
        ligation_site: "Ligation site sequence used for reads trimming."
        max_iter: {
            description: "Maxium number of iterations for the ICE normalization",
            hicpro_field: "MAX_ITER",
        }
        bin_sizes: {
            description: "Resolution of contact maps to generate ",
            hicpro_field: "BIN_SIZE",
        }
        bin_step: {
            description: "Binning step size in `n` coverage i.e. window step."
        }
        matrix_format: {
            description: "Format of the output matrix",
            choices: [
                "complete",
                "asis",
                "upper",
                "lower"
            ],
            hicpro_field: "MATRIX_FORMAT",
        }
        precision: {
            description: "The relative increment in the results before declaring convergence",
            hicpro_field: "EPS",
        }
        filter_low_counts_percentage: {
            description: "Define which percentage of bins with low counts should be force to zero",
            hicpro_field: "FILTER_LOW_COUNT_PERC",
        }
        filter_high_counts_percentage: {
            description: "Define which percentage of bins with low counts should be discarded before normalization",
            hicpro_field: "FILTER_HIGH_COUNT_PERC",
        }
        remove_singleton: {
            description: "Remove singleton reads",
            hicpro_field: "RM_SINGLETON",
        }
        remove_multimapper: {
            description: "Remove multi-mapped reads",
            hicpro_field: "RM_MULTI",
        }
        remove_duplicates: {
            description: "Remove duplicated read pairs",
            hicpro_field: "RM_DUP",
        }
        chromsizes: {
            description: "Tab delimited file with chromosome sizes"
        }
        min_mapq: {
            description: "Minimum mapping quality. Reads with lower quality are discarded.",
            hicpro_field: "MIN_MAPQ",
        }
    }

    input {
        File bowtie_db_tar_gz
        File chromsizes
        Array[File] read_one_fastqs_gz
        String genome_prefix
        String prefix
        File? exclude_list
        File? fragment_file
        File? capture_bed
        File? allele_specific_snp
        String? read_groups
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
        String matrix_format = "upper"
        Boolean remove_singleton = true
        Boolean remove_multimapper = true
        Boolean remove_duplicates = true
        Float filter_low_counts_percentage = 0.02
        Float filter_high_counts_percentage = 0
        Float precision = 0.1
        Int min_mapq = 10
        Int max_iter = 100
        Int bin_step = 1
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

    # scatter (pair in zip(read_one_fastqs_gz, read_two_fastqs_gz)) {
    scatter (pair in zip(flatten(r1_split.fastqs), flatten(r2_split.fastqs))) {
        String fq_prefix = sub(basename(pair.left), ".fastq.gz|.fq.gz|.fastq|.fq", "")
        ReadGroup rg_e2e = ReadGroup{
            ID: "BMG",
            SM: prefix,
        }
        ReadGroup rg_local = ReadGroup{
            ID: "BML",
            SM: prefix,
        }

        # do end-to-end bowtie alignment
        # retain unmapped reads 
        # align read 1, save unaligned reads (--un-gz)
        call bowtie2.align as read1_align { input:
            bowtie_db_tar_gz,
            read_one_fastq_gz = pair.left,
            prefix = fq_prefix + ".R1_global",
            rg = rg_e2e,
            write_unpaired_unaligned = true,
            seed_substring = 30,  # -L
            score_min = "L,-0.6,-0.2",
            end_to_end = true,
            reorder = true,
            max_failed_extends = 20,  # -D
            repetitive_seeds = 3,  # -R
            seed_mismatch = 0,  # -N
            interval_seed_substrings = "S,1,0.50",  # -i
            metrics_file = true,
        }
        # align read 2
        call bowtie2.align as read2_align { input:
            bowtie_db_tar_gz,
            read_one_fastq_gz = pair.right,
            prefix = fq_prefix + ".R2_global",
            rg = rg_e2e,
            write_unpaired_unaligned = true,
            seed_substring = 30,  # -L
            score_min = "L,-0.6,-0.2",
            end_to_end = true,
            reorder = true,
            max_failed_extends = 20,  # -D
            repetitive_seeds = 3,  # -R
            seed_mismatch = 0,  # -N
            interval_seed_substrings = "S,1,0.50",  # -i
            metrics_file = true,
        }

        # ligation site is optional. If not specified, skip local and just use the global alignment.
        if (defined(ligation_site)) {
            FlagFilter unmapped_reads = FlagFilter {
                include_if_all: "0x0",
                exclude_if_any: "0x4",
                include_if_any: "0x0",
                exclude_if_all: "0x0",
            }

            call samtools.filter as read1_filter { input:
                bam = read1_align.aligned_bam,
                bitwise_filter = unmapped_reads,
                prefix = fq_prefix + ".R1_global.mapped",
            }

            call samtools.filter as read2_filter { input:
                bam = read2_align.aligned_bam,
                bitwise_filter = unmapped_reads,
                prefix = fq_prefix + ".R2_global.mapped",
            }

            String site = select_first([ligation_site, ""])
            # Trim FASTQs
            call hicpro.cutsite_trimming as trim_read1 { input:
                fastq = select_first([read1_align.unpaired_unaligned, ""]),
                cutsite = site,
            }
            call hicpro.cutsite_trimming as trim_read2 { input:
                fastq = select_first([read2_align.unpaired_unaligned, ""]),
                cutsite = site,
            }

            # do local bowtie2 alignment
            # align read 1, do not save unaligned reads (no --un-gz)
            call bowtie2.align as read1_local_align { input:
                bowtie_db_tar_gz,
                read_one_fastq_gz = trim_read1.cutsite_trimmed_fastq,
                prefix = fq_prefix + ".R1_local",
                rg = rg_local,
                seed_substring = 20,  # -L
                score_min = "L,-0.6,-0.2",
                end_to_end = true,
                reorder = true,
                max_failed_extends = 20,  # -D
                repetitive_seeds = 3,  # -R
                seed_mismatch = 0,  # -N
                interval_seed_substrings = "S,1,0.50",  # -i
                metrics_file = true,
            }
            # align read 2
            call bowtie2.align as read2_local_align { input:
                bowtie_db_tar_gz,
                read_one_fastq_gz = trim_read2.cutsite_trimmed_fastq,
                prefix = fq_prefix + ".R2_local",
                rg = rg_local,
                seed_substring = 20,  # -L
                score_min = "L,-0.6,-0.2",
                end_to_end = true,
                reorder = true,
                max_failed_extends = 20,  # -D
                repetitive_seeds = 3,  # -R
                seed_mismatch = 0,  # -N
                interval_seed_substrings = "S,1,0.50",  # -i
                metrics_file = true,
            }

            # merge the two alignments
            call samtools.merge as read1_merge { input:
                bams = [read1_filter.filtered_bam, read1_local_align.aligned_bam],
                prefix = fq_prefix + ".R1_merged",
                name_sorted = true,
            }

            call samtools.merge as read2_merge { input:
                bams = [read2_filter.filtered_bam, read2_local_align.aligned_bam],
                prefix = fq_prefix + ".R2_merged",
                name_sorted = true,
            }
        }

        # Sort the merged files
        call samtools.sort as read1_sort { input:
            bam = select_first([read1_merge.merged_bam, read1_align.aligned_bam]),
            natural_name_sort = true,
        }
        call samtools.sort as read2_sort { input:
            bam = select_first([read2_merge.merged_bam, read2_align.aligned_bam]),
            natural_name_sort = true,
        }

        # compute stats
        # total reads, mapped reads, globally mapped reads, locally mapped reads
        call hicpro.mapping_stats as read1_stats { input:
            combined_bam = read1_sort.sorted_bam,
            global_bam = read1_align.aligned_bam,
            local_bam = read1_local_align.aligned_bam,
            prefix = fq_prefix + ".R1",
            tag = "R1",
        }
        call hicpro.mapping_stats as read2_stats { input:
            combined_bam = read2_sort.sorted_bam,
            global_bam = read2_align.aligned_bam,
            local_bam = read2_local_align.aligned_bam,
            prefix = fq_prefix + ".R2",
            tag = "R2",
        }

        call hicpro.bowtie_pairing { input:
            read1_bam = read1_sort.sorted_bam,
            read2_bam = read2_sort.sorted_bam,
            prefix = fq_prefix,
            remove_singleton,
            remove_multimapper,
        }

        call hicpro.mapped_2hic_fragments { input:
            mapped_reads = bowtie_pairing.combined_bam,
            fragment = fragment_file,
            addl_output = true,
            sam = true,
        }
    }

    call hicpro.merge_valid_interactions { input:
        interactions = mapped_2hic_fragments.valid_pairs,
        prefix,
        remove_duplicates,
    }

    call hicpro.merge_stats { input:
        read1_mapping_stats = read1_stats.mapping_stats,
        read2_mapping_stats = read2_stats.mapping_stats,
        valid_pairs_stats = bowtie_pairing.combined_stats,
        rs_stats = mapped_2hic_fragments.rs_stats,
        prefix,
    }

    call hicpro.build_raw_maps { input:
        hic_file = merge_valid_interactions.all_valid_pairs,
        bin_sizes,
        genome_fragment = fragment_file,
        matrix_format,
        chromsizes,
        prefix,
    }

    call hicpro.qc_hic { input:
        plot_type = "all",
        mapping_stats = flatten([
            read1_stats.mapping_stats,
            read2_stats.mapping_stats,
        ]),
        pairing_stats = bowtie_pairing.combined_stats,
        fragment_stats = flatten([mapped_2hic_fragments.rs_stats, mapped_2hic_fragments.valid_pairs]),
        contacts_stats = [merge_valid_interactions.all_valid_pairs_stats],
        sample_name = prefix,
        remove_singleton,
        remove_multimapper,
    }

    call hicpro.ice_normalization { input:
        bin_sizes,
        contact_counts = build_raw_maps.contact_counts,
        prefix,
    }

    call hilow.converthic { input:
        all_valid_pairs = merge_valid_interactions.all_valid_pairs,
        chromsizes,
    }

    if (defined(exclude_list)) {
        call hilow.filter { input:
            all_valid_pairs = merge_valid_interactions.all_valid_pairs,
            chromsizes,
            exclude_list = select_first([exclude_list, ""]),
        }
    }

    call hilow.qcreport { input:
        prefix,
        all_valid_pairs_stats = merge_valid_interactions.all_valid_pairs_stats,
        mapping_stats_R1 = merge_stats.read1_mapping_stats_merged,
        mapping_stats_R2 = merge_stats.read2_mapping_stats_merged,
        pairing_stats = merge_stats.valid_pairs_stats_merged,
    }

    output {
        File hic_file = converthic.hic_file
        File? filtered_pairs = filter.filtered_pairs
        File all_valid_pairs = merge_valid_interactions.all_valid_pairs
        File qc_report = qcreport.qc_report
        File? mapping_stats_plot = qc_hic.mapping_stats_plot
        File? pairing_stats_plot = qc_hic.pairing_stats_plot
        File? filtering_stats_plot = qc_hic.filtering_stats_plot
        File? filtering_size_plot = qc_hic.filtering_size_plot
        File? contacts_stats_plot = qc_hic.contacts_stats_plot
        Array[File] ice_normalized_matrices = ice_normalization.iced_matrices
        Array[File] combined_bams = bowtie_pairing.combined_bam
    }
}
