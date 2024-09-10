version 1.1

import "../../data_structures/flag_filter.wdl"
import "../../data_structures/read_group.wdl"
import "../../tools/bowtie2.wdl"
import "../../tools/hilow.wdl"
import "../../tools/samtools.wdl"
import "../../tools/util.wdl"

workflow hicpro_core {
    meta {
        description: "HiC-Pro implementation."
        outputs: {
            all_valid_pairs: "All valid pairs file",
            ice_normalized_matrices: "ICE normalized matrices",
            combined_bams: "Combined BAM files with gloabl and (if run) local alignments",
            mapping_stats_plot: "Mapping stats plot for R1 and R2 tags. Contains aligned percentage, unaligned percentage, and full-read and trimmed read mapping.",
            pairing_stats_plot: "Pairing stats plot containing reported status for all pairs and quality for filtered pairs.",
            filtering_stats_plot: "Plot of pair alignments to restriction fragments. Contains orientation and classification of intercations.",
            filtering_size_plot: "Plot of distribution of fragment sizes.",
            contacts_stats_plot: "Plot of contact ranges for valid pairs.",
            contact_stats: "Contact statistics",
            merged_read_one_mapping_stats: "Merged read 1 mapping statistics",
            merged_read_two_mapping_stats: "Merged read 2 mapping statistics",
            merged_pairing_stats: "Merged pairing statistics",
        }
    }

    parameter_meta {
        read_one_fastqs_gz: "An array of gzipped FASTQ files containing read one information"
        read_groups: "An array of ReadGroup structs containing read group information for each input FASTQ to output in the BAM file"
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
        Array[ReadGroup] read_groups
        String prefix
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

    scatter (tuple in zip(zip(read_one_fastqs_gz, read_two_fastqs_gz), read_groups)) {
        call util.split_fastq as r1_split { input:
            fastq = tuple.left.left,
            reads_per_file = 30000000,
        }
        call util.split_fastq as r2_split { input:
            fastq = tuple.left.right,
            reads_per_file = 30000000,
        }

        scatter (pair in zip(r1_split.fastqs, r2_split.fastqs)) {
            String fq_prefix = sub(basename(pair.left), ".fastq.gz|.fq.gz|.fastq|.fq", "")
            # do end-to-end bowtie alignment
            # retain unmapped reads
            # align read 1, save unaligned reads (--un-gz)
            call bowtie2.align as read1_align { input:
                bowtie_db_tar_gz,
                read_one_fastq_gz = pair.left,
                prefix = fq_prefix + ".R1_global",
                rg = tuple.right,
                write_unpaired_unaligned = true,
                seed_substring = 30,  # -L
                score_min = {
                    "function_type": "L",
                    "constant": -0.6,
                    "coefficient": -0.2,
                },
                end_to_end = true,
                reorder = true,
                max_failed_extends = 20,  # -D
                repetitive_seeds = 3,  # -R
                seed_mismatch = 0,  # -N
                interval_seed_substrings = {
                    "function_type": "S",
                    "constant": 1,
                    "coefficient": 0.50,
                },  # -i
                metrics_file = true,
            }
            # align read 2
            call bowtie2.align as read2_align { input:
                bowtie_db_tar_gz,
                read_one_fastq_gz = pair.right,
                prefix = fq_prefix + ".R2_global",
                rg = tuple.right,
                write_unpaired_unaligned = true,
                seed_substring = 30,  # -L
                score_min = {
                    "function_type": "L",
                    "constant": -0.6,
                    "coefficient": -0.2,
                },
                end_to_end = true,
                reorder = true,
                max_failed_extends = 20,  # -D
                repetitive_seeds = 3,  # -R
                seed_mismatch = 0,  # -N
                interval_seed_substrings = {
                    "function_type": "S",
                    "constant": 1,
                    "coefficient": 0.50,
                },  # -i
                metrics_file = true,
            }

            # ligation site is optional.
            # If not specified, skip local and just use the global alignment.
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
                call cutsite_trimming as trim_read1 { input:
                    fastq = select_first([read1_align.unpaired_unaligned, ""]),
                    cutsite = site,
                }
                call cutsite_trimming as trim_read2 { input:
                    fastq = select_first([read2_align.unpaired_unaligned, ""]),
                    cutsite = site,
                }

                # do local bowtie2 alignment
                # align read 1, do not save unaligned reads (no --un-gz)
                call bowtie2.align as read1_local_align { input:
                    bowtie_db_tar_gz,
                    read_one_fastq_gz = trim_read1.cutsite_trimmed_fastq,
                    prefix = fq_prefix + ".R1_local",
                    rg = tuple.right,
                    seed_substring = 20,  # -L
                    score_min = {
                        "function_type": "L",
                        "constant": -0.6,
                        "coefficient": -0.2,
                    },
                    end_to_end = true,
                    reorder = true,
                    max_failed_extends = 20,  # -D
                    repetitive_seeds = 3,  # -R
                    seed_mismatch = 0,  # -N
                    interval_seed_substrings = {
                        "function_type": "S",
                        "constant": 1,
                        "coefficient": 0.50,
                    },  # -i
                    metrics_file = true,
                }
                # align read 2
                call bowtie2.align as read2_local_align { input:
                    bowtie_db_tar_gz,
                    read_one_fastq_gz = trim_read2.cutsite_trimmed_fastq,
                    prefix = fq_prefix + ".R2_local",
                    rg = tuple.right,
                    seed_substring = 20,  # -L
                    score_min = {
                        "function_type": "L",
                        "constant": -0.6,
                        "coefficient": -0.2,
                    },
                    end_to_end = true,
                    reorder = true,
                    max_failed_extends = 20,  # -D
                    repetitive_seeds = 3,  # -R
                    seed_mismatch = 0,  # -N
                    interval_seed_substrings = {
                        "function_type": "S",
                        "constant": 1,
                        "coefficient": 0.50,
                    },  # -i
                    metrics_file = true,
                }

                # merge the two alignments
                call samtools.merge as read1_merge { input:
                    bams = [read1_filter.filtered_bam, read1_local_align.aligned_bam],
                    prefix = fq_prefix + ".R1_merged",
                    name_sorted = true,
                    attach_rg = false,
                }

                call samtools.merge as read2_merge { input:
                    bams = [read2_filter.filtered_bam, read2_local_align.aligned_bam],
                    prefix = fq_prefix + ".R2_merged",
                    name_sorted = true,
                    attach_rg = false,
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
            call mapping_stats as read1_stats { input:
                combined_bam = read1_sort.sorted_bam,
                global_bam = read1_align.aligned_bam,
                local_bam = read1_local_align.aligned_bam,
                prefix = fq_prefix + ".R1",
                tag = "R1",
            }
            call mapping_stats as read2_stats { input:
                combined_bam = read2_sort.sorted_bam,
                global_bam = read2_align.aligned_bam,
                local_bam = read2_local_align.aligned_bam,
                prefix = fq_prefix + ".R2",
                tag = "R2",
            }

            call bowtie_pairing { input:
                read1_bam = read1_sort.sorted_bam,
                read2_bam = read2_sort.sorted_bam,
                prefix = fq_prefix,
                remove_singleton,
                remove_multimapper,
            }

            call mapped_2hic_fragments { input:
                mapped_reads = bowtie_pairing.combined_bam,
                fragment = fragment_file,
                addl_output = true,
                sam = true,
            }
        }
    }

    call merge_valid_interactions { input:
        interactions = flatten(mapped_2hic_fragments.valid_pairs),
        prefix,
        remove_duplicates,
    }

    call merge_stats { input:
        read1_mapping_stats = flatten(read1_stats.mapping_stats),
        read2_mapping_stats = flatten(read2_stats.mapping_stats),
        valid_pairs_stats = flatten(bowtie_pairing.combined_stats),
        rs_stats = flatten(mapped_2hic_fragments.rs_stats),
        prefix,
    }

    call qc_hic { input:
        plot_type = "all",
        mapping_stats = flatten([
            flatten(read1_stats.mapping_stats),
            flatten(read2_stats.mapping_stats),
        ]),
        pairing_stats = flatten(bowtie_pairing.combined_stats),
        fragment_stats = flatten([
            flatten(mapped_2hic_fragments.rs_stats),
            flatten(mapped_2hic_fragments.valid_pairs),
        ]),
        contacts_stats = [merge_valid_interactions.all_valid_pairs_stats],
        sample_name = prefix,
        remove_singleton = true,
        remove_multimapper = true,
    }

    call build_raw_maps { input:
        hic_file = merge_valid_interactions.all_valid_pairs,
        bin_sizes,
        genome_fragment = fragment_file,
        matrix_format,
        chromsizes,
        prefix,
    }

    call ice_normalization { input:
        bin_sizes,
        contact_counts = build_raw_maps.contact_counts,
        prefix,
    }

    output {
        File all_valid_pairs = merge_valid_interactions.all_valid_pairs
        Array[File] ice_normalized_matrices = ice_normalization.iced_matrices
        Array[File] combined_bams = flatten(bowtie_pairing.combined_bam)
        File? mapping_stats_plot = qc_hic.mapping_stats_plot
        File? pairing_stats_plot = qc_hic.pairing_stats_plot
        File? filtering_stats_plot = qc_hic.filtering_stats_plot
        File? filtering_size_plot = qc_hic.filtering_size_plot
        File? contacts_stats_plot = qc_hic.contacts_stats_plot
        File contact_stats = merge_valid_interactions.all_valid_pairs_stats
        File merged_read_one_mapping_stats = merge_stats.read1_mapping_stats_merged
        File merged_read_two_mapping_stats = merge_stats.read2_mapping_stats_merged
        File merged_pairing_stats = merge_stats.valid_pairs_stats_merged
    }
}

task cutsite_trimming {
    meta {
        description: "Trim cutsite from reads"
        outputs: {
            cutsite_trimmed_fastq: "FASTQ file with cutsites trimmed",
            cutsite_trimmed_log: "Log file of the trimming process",
        }
    }

    parameter_meta {
        fastq: "Input FASTQ file"
        cutsite: "Cutsite to trim"
    }

    input {
        File fastq
        String cutsite
    }

    String prefix = sub(basename(fastq), ".fq.gz|.fastq.gz|.fq|.fastq", "")

    command <<<
        set -euo pipefail

        fastq=~{basename(fastq, ".gz")}
        gunzip -c ~{fastq} > "$fastq" \
            || ln -sf ~{fastq} "$fastq"

        /HiC-Pro_3.0.0/scripts/cutsite_trimming \
            --fastq $fastq \
            --cutsite ~{cutsite} \
            --out ~{prefix}_cutsite_trimmed.fastq \
            > ~{prefix}_readsTrimming.log \
            2>&1

        gzip ~{prefix}_cutsite_trimmed.fastq && rm $fastq
    >>>

    output {
        File cutsite_trimmed_fastq = prefix + "_cutsite_trimmed.fastq.gz"
        File cutsite_trimmed_log = prefix + "_readsTrimming.log"
    }

    runtime {
        container: "nservant/hicpro:3.0.0"
        cpu: 1
        memory: "4 GB"
        maxRetries: 1
    }
}

task bowtie_pairing {
    meta {
        description: "Merge R1 and R2 BAM files into one paired-end BAM file"
        outputs: {
            combined_bam: "Combined BAM file",
            allele_specific_bam: "Allele specific BAM file",
            combined_stats: "Paired-end statistics file",
        }
    }

    parameter_meta {
        read1_bam: "BAM file containing read 1"
        read2_bam: "BAM file containing read 2"
        allele_specific_tag: "Tag to add to the output"
        prefix: "Prefix for the output file"
        remove_singleton: "Remove singleton reads"
        remove_multimapper: "Remove multi-mapped reads"
        min_mapq: "Minimum mapping quality"
    }

    input {
        File read1_bam
        File read2_bam
        File? allele_specific_tag
        String prefix = basename(read1_bam, ".bwt2merged.bam")
        Boolean remove_singleton = false
        Boolean remove_multimapper = false
        Int min_mapq = 10
    }

    command <<<
        set -euo pipefail
        # merge pairs
        python /HiC-Pro_3.0.0/scripts/mergeSAM.py \
            -q ~{min_mapq} \
            -t \
            -v \
            ~{if remove_singleton then "" else "-s"} \
            ~{if remove_multimapper then "" else "-m"} \
            -f ~{read1_bam} \
            -r ~{read2_bam} \
            -o ~{prefix}.bwt2pairs.bam

        ~{(
            if defined(allele_specific_tag)
            then "python /HiC-Pro_3.0.0/scripts/addTagToBAM.py -s ~{allele_specific_tag} -v -r -i ~{prefix}.bwt2pairs.bam -o ~{prefix}.bwt2pairs_allspe.bam"
            else ""
        )}

    >>>

    output {
        File combined_bam = prefix + ".bwt2pairs.bam"
        File? allele_specific_bam = prefix + ".bwt2pairs_allspe.bam"
        File combined_stats = prefix + ".bwt2pairs.pairstat"
    }

    runtime {
        container: "nservant/hicpro:3.0.0"
        maxRetries: 1
    }
}

task mapping_stats {
    meta {
        description: "Compute mapping statistics"
        outputs: {
            mapping_stats: "Summarized mapping statistics file"
        }
    }

    parameter_meta {
        combined_bam: "Combined BAM file"
        global_bam: "End-to-end alignments BAM file"
        local_bam: "Local alignments BAM file"
        tag: "Tag to add to the output"
        prefix: "Prefix for the output file"
    }

    input {
        File combined_bam
        File global_bam
        File? local_bam
        String? tag
        String prefix = basename(combined_bam, ".bwt2glob.bam")
    }

    command <<<
        set -euo pipefail

        echo "## HiC-Pro Mapping Statistics" > ~{prefix}.mapstat
        echo "## ~{prefix}.mapstat" >> ~{prefix}.mapstat

        total_reads=$(samtools view -c ~{combined_bam})

        mapped_reads=$(samtools view -c -F 4 ~{combined_bam})

        global_mapped_reads=$(samtools view -c -F 4 ~{global_bam})

        if ~{if defined(local_bam) then true else false}
        then
            local_mapped_reads=$(samtools view -c -F 4 ~{local_bam})
        fi

        echo -e "total_~{tag}\t$total_reads" >> ~{prefix}.mapstat
        echo -e "mapped_~{tag}\t$mapped_reads" >> ~{prefix}.mapstat
        echo -e "global_~{tag}\t$global_mapped_reads" >> ~{prefix}.mapstat
        echo -e "local_~{tag}\t$local_mapped_reads" >> ~{prefix}.mapstat
    >>>

    output {
        File mapping_stats = prefix + ".mapstat"
    }

    runtime {
        container: "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"
        maxRetries: 1
    }
}

task mapped_2hic_fragments {
    meta {
        description: "Keep only valid 3C products"
        outputs: {
            valid_pairs: "Valid pairs file",
            rs_stats: "RS statistics file",
        }
    }

    parameter_meta {
        fragment: "Restriction fragment file in GFF3 format"
        mapped_reads: "Mapped reads in BAM/SAM format"
        shortest_insert_size: {
            description: "Shortest insert size of mapped reads to consider",
            hicpro_field: "MIN_INSERT_SIZE",
        }
        longest_insert_size: {
            description: "Longest insert size of mapped reads to consider",
            hicpro_field: "MAX_INSERT_SIZE",
        }
        shortest_fragment_length: {
            description: "Shortest restriction fragment length to consider",
            hicpro_field: "MIN_FRAG_SIZE",
        }
        longest_fragment_length: {
            description: "Longest restriction fragment length to consider",
            hicpro_field: "MAX_FRAG_SIZE",
        }
        min_cis_distance: {
            description: "Minimum distance between intrachromosomal contact to consider. Filters contacts below this distance. Mainly useful for DNase Hi-C",
            hicpro_field: "MIN_CIS_DIST",
        }
        genotype_tag: "Adds XA tag to report in the valid pairs output for allele specific classification"
        addl_output: {
            description: "Write all additional output files, with information about the discarded reads (self-circle, dangling end, etc.)",
            hicpro_field: "GET_ALL_INTERACTION_CLASSES",
        }
        sam: {
            description: "Output an additional SAM file with flag 'CT' for pairs classification",
            hicpro_field: "GET_PROCESS_BAM",
        }
    }

    input {
        File mapped_reads
        File? fragment
        Boolean genotype_tag = false
        Boolean addl_output = false
        Boolean sam = false
        Int shortest_insert_size = 0
        Int longest_insert_size = 0
        Int shortest_fragment_length = 0
        Int longest_fragment_length = 0
        Int min_cis_distance = 0
    }

    String outfile_name = basename(mapped_reads, ".bam") + ".validPairs"

    command <<<
        set -euo pipefail

        if ~{if defined(fragment) then true else false}
        then
            frag=~{basename(select_first([fragment, ""]), ".gz")}
            gunzip -c ~{fragment} > "$frag" \
                || ln -sf ~{fragment} "$frag"
        fi

        ~{(
            if defined(fragment)
            then "python /HiC-Pro_3.0.0/scripts/mapped_2hic_fragments.py -f $frag"  # Fragment file found
            else "python /HiC-Pro_3.0.0/scripts/mapped_2hic_dnase.py"  # DNAse
        )} \
        -v \
        -r ~{mapped_reads} \
        ~{if defined(fragment) && sam then "-S" else ""} \
        ~{if addl_output then "-a" else ""} \
        ~{if min_cis_distance > 0 then "-d " + min_cis_distance else ""} \
        ~{if genotype_tag then "-g XA" else ""} \
        ~{(
            if defined(fragment) && shortest_insert_size > 0
            then "-s " + shortest_insert_size
            else ""
        )} \
        ~{(
            if defined(fragment) && longest_insert_size > 0
            then "-l " + longest_insert_size
            else ""
        )} \
        ~{(
            if defined(fragment) && shortest_fragment_length > 0
            then "-t " + shortest_fragment_length
            else ""
        )} \
        ~{(
            if defined(fragment) && longest_fragment_length > 0
            then "-m " + longest_fragment_length
            else ""
        )}

        LANG=en; sort -k2,2V -k3,3n -k5,5V -k6,6n -o ~{outfile_name} ~{outfile_name}
    >>>

    output {
        File valid_pairs = outfile_name
        File rs_stats = basename(outfile_name, ".validPairs") + ".RSstat"
    }

    runtime {
        container: "nservant/hicpro:3.0.0"
        maxRetries: 1
    }
}

task merge_valid_interactions {
    meta {
        description: "Merge valid pairs files"
        outputs: {
            all_valid_pairs: "Merged valid pairs file",
            on_target: "On target valid pairs file",
            split_interactions: "Split interactions file",
            all_valid_pairs_stats: "Merged valid pairs statistics file",
        }
    }

    parameter_meta {
        interactions: "Valid pairs files to merge"
        prefix: "Prefix for the output files"
        remove_duplicates: "Remove duplicates"
        report_capture: "Report capture"
        allele_specific_snp: "Use allele specific SNP"
        capture_target: "Capture target file"
    }

    input {
        Array[File] interactions
        String prefix
        String? capture_target
        Boolean remove_duplicates = false
        Boolean report_capture = false
        Boolean allele_specific_snp = false
    }

    command <<<
        set -euo pipefail

        if ~{remove_duplicates}
        then
            LANG=en; sort \
                -k2,2V -k3,3n -k5,5V -k6,6n \
                -m ~{sep(" ", interactions)} \
                | awk -F "\t" \
                'BEGIN{c1=0;c2=0;s1=0;s2=0}(c1!=$2 || c2!=$5 || s1!=$3 || s2!=$6){print;c1=$2;c2=$5;s1=$3;s2=$6}' \
                > ~{prefix}.allValidPairs
        else
            cat ~{sep(" ", interactions)} > ~{prefix}.allValidPairs
        fi

        allcount=$(cat ~{sep(" ", interactions)} | wc -l)
        allcount_rmdup=$(cat ~{prefix}.allValidPairs | wc -l)
        ndbup=$(( $allcount - $allcount_rmdup ))

        ## merge stat file
        echo -e "valid_interaction\t"$allcount > ~{prefix}_allValidPairs.mergestat
        echo -e "valid_interaction_rmdup\t"$allcount_rmdup >> ~{prefix}_allValidPairs.mergestat
        awk 'BEGIN{cis=0;trans=0;sr=0;lr=0} $2 == $5{cis=cis+1; d=$6>$3?$6-$3:$3-$6; if (d<=20000){sr=sr+1}else{lr=lr+1}} $2!=$5{trans=trans+1}END{print "trans_interaction\t"trans"\ncis_interaction\t"cis"\ncis_shortRange\t"sr"\ncis_longRange\t"lr}' ~{prefix}.allValidPairs >> ~{prefix}_allValidPairs.mergestat

        if ~{if defined(capture_target) then "true" else "false"}
        then
            python /HiC-Pro_3.0.0/scripts/onTarget.py \
                -i ~{prefix}.allValidPairs \
                -t ~{capture_target} \
                ~{if (report_capture) then "--cis" else ""} \
                -s ~{prefix}_allValidPairs.mergestat \
                -v \
                > ~{prefix}_ontarget.allValidPairs
        fi

        if ~{allele_specific_snp}
        then
            python /HiC-Pro_3.0.0/scripts/split_valid_interactions.py \
                ~{(
                    if defined(capture_target)
                    then "-i ~{prefix}_ontarget.allValidPairs"
                    else "-i ~{prefix}.allValidPairs"
                )} \
                -s ~{prefix}_allValidPairs_assplit.stat \
                -v    
        fi
    >>>

    output {
        File all_valid_pairs = prefix + ".allValidPairs"
        File? on_target = prefix + "_ontarget.allValidPairs"
        File? split_interactions = prefix + "_allValidPairs_assplit.stat"
        File all_valid_pairs_stats = prefix + "_allValidPairs.mergestat"
    }

    runtime {
        container: "nservant/hicpro:3.0.0"
        maxRetries: 1
    }
}

task merge_stats {
    meta {
        description: "Merge statistics files"
        outputs: {
            read1_mapping_stats_merged: "Mapping statistics for read 1",
            read2_mapping_stats_merged: "Mapping statistics for read 2",
            valid_pairs_stats_merged: "Valid pairs statistics",
            rs_stats_merged: "RS statistics",
        }
    }

    parameter_meta {
        read1_mapping_stats: "Mapping statistics for read 1"
        read2_mapping_stats: "Mapping statistics for read 2"
        valid_pairs_stats: "Valid pairs statistics"
        rs_stats: "RS statistics"
        prefix: "Prefix for the output files"
    }

    input {
        Array[File] read1_mapping_stats
        Array[File] read2_mapping_stats
        Array[File] valid_pairs_stats
        Array[File] rs_stats
        String prefix
    }

    command <<<
        set -euo pipefail

        # merge read1 mapping stats
        mkdir read1
        ln -s ~{sep(" ", read1_mapping_stats)} read1/
        python /HiC-Pro_3.0.0/scripts/merge_statfiles.py \
            -d read1 \
            -p "*.R1*.mapstat" \
            -v \
            > ~{prefix}_R1.mmapstat
        # merge read2 mapping stats
        mkdir read2
        ln -s ~{sep(" ", read2_mapping_stats)} read2/
        python /HiC-Pro_3.0.0/scripts/merge_statfiles.py \
            -d read2 \
            -p "*.R2*.mapstat" \
            -v \
            > ~{prefix}_R2.mmapstat
        # merge pairing stats
        mkdir pairs
        ln -s ~{sep(" ", valid_pairs_stats)} pairs/
        python /HiC-Pro_3.0.0/scripts/merge_statfiles.py \
            -d pairs \
            -p "*.pairstat" \
            -v \
            > ~{prefix}.mpairstat
        # merge RS stat
        mkdir rsstat
        ln -s ~{sep(" ", rs_stats)} rsstat/
        python /HiC-Pro_3.0.0/scripts/merge_statfiles.py \
            -d rsstat \
            -p "*.RSstat" \
            -v \
            > ~{prefix}.mRSstat
    >>>

    output {
        File read1_mapping_stats_merged = prefix + "_R1.mmapstat"
        File read2_mapping_stats_merged = prefix + "_R2.mmapstat"
        File valid_pairs_stats_merged = prefix + ".mpairstat"
        File rs_stats_merged = prefix + ".mRSstat"
    }

    runtime {
        container: "nservant/hicpro:3.0.0"
        maxRetries: 1
    }
}

task build_raw_maps {
    meta {
        description: "Build raw maps from valid pairs for all specified resolutions"
        outputs: {
            contact_counts: "Contact counts files for each resolution",
            contact_counts_bed: "Contact counts for each resolution in BED format",
        }
    }

    parameter_meta {
        chromsizes: "Chromosome sizes file"
        hic_file: "Hi-C file"
        bin_sizes: "Bin sizes for each contact counts file"
        prefix: "Prefix for the output files"
        allele_specific_snp: "Use allele specific SNP"
        capture_target: "Capture target file"
        genome_fragment: "Genome fragment file"
        matrix_format: "Matrix format"
    }

    input {
        File chromsizes
        File hic_file
        Array[Int] bin_sizes
        String prefix
        File? genome_fragment
        String matrix_format = "upper"
        Boolean allele_specific_snp = false
        Boolean capture_target = false
    }

    command <<<
        set -euo pipefail

        for bin_size in ~{sep(" ", bin_sizes)}
        do
            cat ~{hic_file} \
            | /HiC-Pro_3.0.0/scripts/build_matrix \
                --matrix-format ~{matrix_format} \
                ~{(
                    if (length(bin_sizes) > 1)
                    then "--binsize $bin_size"
                    else "--binfile ~{select_first([genome_fragment, ""])}"
                )} \
                --chrsizes ~{chromsizes} \
                --ifile /dev/stdin \
                --oprefix ~{prefix}_${bin_size}
        done
    >>>

    output {
        Array[File] contact_counts = glob(prefix + "_*.matrix")
        Array[File] contact_counts_bed = glob(prefix + "_*.bed")
    }

    runtime {
        container: "nservant/hicpro:3.0.0"
        maxRetries: 1
    }
}

task ice_normalization {
    meta {
        description: "Normalize Hi-C data using the ICE algorithm to correct several sources of bias"
        external_help: "https://nservant.github.io/HiC-Pro/MANUAL.html"
        outputs: {
            iced_matrices: "ICE normalized matrices",
            biases: "Biases files",
        }
    }

    parameter_meta {
        contact_counts: "Contact counts files to normalize"
        bin_sizes: "Bin sizes for each contact counts file"
        prefix: "Prefix for the output files"
        filtering_percentage: "Percentage of contacts to filter out"
        remove_all_zeroes_loci: "Remove all zero loci"
        filter_low_counts_percentage: "Filter low counts percentage"
        filter_high_counts_percentage: "Filter high counts percentage"
        precision: "Precision of the normalization"
        max_iterations: "Maximum number of iterations"
    }

    input {
        Array[File] contact_counts
        Array[Int] bin_sizes
        String prefix
        Float? filtering_percentage
        Boolean remove_all_zeroes_loci = false
        Float filter_low_counts_percentage = 0.02
        Float filter_high_counts_percentage = 0.0
        Float precision = 0.1
        Int max_iterations = 100
    }

    command <<<
        for map in ~{sep(" ", contact_counts)}
        do
            name=$(basename $map)
            for bin in ~{sep(" ", bin_sizes)}
            do
                ice \
                    --results_filename ${name}_${bin}_iced.matrix \
                    --filter_low_counts_perc ~{filter_low_counts_percentage} \
                    --filter_high_counts_perc ~{filter_high_counts_percentage} \
                    --max_iter ~{max_iterations} \
                    --eps ~{precision} \
                    --remove-all-zeros-loci \
                    --output-bias 1 \
                    $map
            done
        done
    >>>

    output {
        Array[File] iced_matrices = glob(prefix + "*_iced.matrix")
        Array[File] biases = glob(prefix + "*_iced.matrix.biases")
    }

    runtime {
        container: "nservant/hicpro:3.0.0"
        maxRetries: 1
        cpu: 1
        memory: "4 GB"
    }
}

task qc_hic {
    meta {
        description: "Plot Hi-C quality control statistics"
        outputs: {
            mapping_stats_plot: "Plot of mapping statistcs",
            pairing_stats_plot: "Plot of fragment pairing",
            filtering_stats_plot: "Plot of Hi-C fragments",
            filtering_size_plot: "Plot of Hi-C fragment sizes",
            contacts_stats_plot: "Plot of Hi-C contact ranges",
        }
    }

    parameter_meta {
        mapping_stats: "Mapping statistics files"
        pairing_stats: "Pairing statistics files"
        fragment_stats: "Fragment statistics files"
        contacts_stats: "Contacts statistics files"
        plot_type: {
            description: "Type of plot to generate",
            choices: [
                "all",
                "mapping",
                "pairing",
                "filtering",
                "contacts"
            ],
        }
        sample_name: "Sample name to use in plot labels"
        remove_singleton: "Remove singleton reads"
        remove_multimapper: "Remove multi-mapped reads"
    }

    input {
        Array[File] mapping_stats
        Array[File] pairing_stats
        Array[File] fragment_stats
        Array[File] contacts_stats
        String plot_type = "all"
        String sample_name = "sample"
        Boolean remove_singleton = true
        Boolean remove_multimapper = true
    }

    Int rm_single_arg = if remove_singleton then 1 else 0
    Int rm_multi_arg = if remove_multimapper then 1 else 0

    command <<<
        set -euo pipefail

        mkdir plots

        mkdir bowtie
        ln -sf ~{sep(" ", mapping_stats)} bowtie/
        ln -sf ~{sep(" ", pairing_stats)} bowtie/

        mkdir hic
        ln -sf ~{sep(" ", fragment_stats)} hic/

        mkdir stats
        ln -sf ~{sep(" ", contacts_stats)} stats/

        # mapping plots
        if [[ "~{plot_type}" == "all" || "~{plot_type}" == "mapping" ]]
        then
            R CMD BATCH \
                --no-save \
                --no-restore \
                "--args picDir='plots' bwtDir='bowtie' sampleName='~{sample_name}' r1tag='.R1' r2tag='.R2'" \
                /HiC-Pro_3.0.0/scripts/plot_mapping_portion.R \
                plot_mapping_portion.Rout
        fi

        # pairing plots
        if [[ "~{plot_type}" == "all" || "~{plot_type}" == "pairing" ]]
        then
            R CMD BATCH \
                --no-save \
                --no-restore \
                "--args picDir='plots' bwtDir='bowtie' sampleName='~{sample_name}' rmMulti='~{rm_multi_arg}' rmSingle='~{rm_single_arg}'" \
                /HiC-Pro_3.0.0/scripts/plot_pairing_portion.R \
                plot_pairing_portion.Rout
        fi
        
        # filtering plots
        if [[ "~{plot_type}" == "all" || "~{plot_type}" == "filtering" ]]
        then
            R CMD BATCH \
                --no-save \
                --no-restore \
                "--args picDir='plots' hicDir='hic' sampleName='~{sample_name}' rmMulti='~{rm_multi_arg}' rmSingle='~{rm_single_arg}'" \
                /HiC-Pro_3.0.0/scripts/plot_hic_fragment.R \
                plot_hic_fragment.Rout
        fi

        # contacts plots
        if [[ "~{plot_type}" == "all" || "~{plot_type}" == "contacts" ]]
        then
            R CMD BATCH \
                --no-save \
                --no-restore \
                "--args picDir='plots' hicDir='hic' statsDir='stats' sampleName='~{sample_name}'" \
                /HiC-Pro_3.0.0/scripts/plot_hic_contacts.R \
                plot_hic_contacts.Rout
        fi
    >>>

    output {
        File? mapping_stats_plot = glob("plots/plotMapping_*.pdf")[0]
        File? pairing_stats_plot = glob("plots/plotMappingPairing_*.pdf")[0]
        File? filtering_stats_plot = glob("plots/plotHiCFragment_*.pdf")[0]
        File? filtering_size_plot = glob("plots/plotHiCFragmentSize_*.pdf")[0]
        File? contacts_stats_plot = glob("plots/plotHiCContactRanges_*.pdf")[0]
    }

    runtime {
        container: "nservant/hicpro:3.0.0"
        cpu: 1
        memory: "4 GB"
        maxRetries: 1
    }
}
