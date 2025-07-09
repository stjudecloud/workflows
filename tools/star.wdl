## [Homepage](https://github.com/alexdobin/STAR)

version 1.1

task build_star_db {
    meta {
        description: "Runs STAR's build command to generate a STAR format reference for alignment"
        outputs: {
            star_db: "A gzipped TAR file containing the STAR reference files. Suitable as the `star_db_tar_gz` input to the `alignment` task."
        }
    }

    parameter_meta {
        reference_fasta: "The FASTA format reference file for the genome"
        gtf: "GTF format feature file"
        db_name: {
            description: "Name for output in compressed, archived format. The suffix `.tar.gz` will be added.",
            group: "Common",
        }
        sjdb_gtf_chr_prefix: {
            description: "prefix for chromosome names in a GTF file (e.g. 'chr' for using ENSEMBL annotations with UCSC genomes)",
            group: "Common",
        }
        sjdb_gtf_feature_exon: "feature type in GTF file to be used as exons for building transcripts"
        sjdb_gtf_tag_exon_parant_transcript: "GTF attribute name for parent transcript ID"
        sjdb_gtf_tag_exon_parent_gene: "GTF attribute name for parent gene ID"
        sjdb_gtf_tag_exon_parent_gene_name: "GTF attribute name for parent gene name"
        sjdb_gtf_tag_exon_parent_gene_type: "GTF attribute name for parent gene type"
        use_all_cores: {
            description: "Use all cores? Recommended for cloud environments.",
            group: "Resources",
        }
        genome_chr_bin_n_bits: "=log2(chrBin), where chrBin is the size of the bins for genome storage: each chromosome will occupy an integer number of bins. For a genome with large number of contigs, it is recommended to scale this parameter as min(18, log2[max(GenomeLength/NumberOfReferences,ReadLength)])."
        genome_SA_index_n_bases: "length (bases) of the SA pre-indexing string. Typically between 10 and 15. Longer strings will use much more memory, but allow faster searches. For small genomes, the parameter `--genomeSAindexNbases` must be scaled down to `min(14, log2(GenomeLength)/2 - 1)`."
        genome_SA_sparse_d: "suffix array sparsity, i.e. distance between indices: use bigger numbers to decrease needed RAM at the cost of mapping speed reduction."
        genome_suffix_length_max: "maximum length of the suffixes, has to be longer than read length. -1 = infinite."
        sjdb_overhang: {
            description: "length of the donor/acceptor sequence on each side of the junctions, ideally = (mate_length - 1). **[STAR default]**: `100`. **[WDL default]**: `125`.",
            group: "Common",
        }
        ncpu: {
            description: "Number of cores to allocate for task",
            group: "Resources",
        }
        memory_gb: {
            description: "RAM to allocate for task, specified in GB",
            group: "Resources",
        }
        modify_disk_size_gb: {
            description: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB.",
            group: "Resources",
        }
    }

    input {
        File reference_fasta
        File gtf
        String db_name = "star_db"
        String sjdb_gtf_chr_prefix = "-"
        String sjdb_gtf_feature_exon = "exon"
        String sjdb_gtf_tag_exon_parant_transcript = "transcript_id"
        String sjdb_gtf_tag_exon_parent_gene = "gene_id"
        String sjdb_gtf_tag_exon_parent_gene_name = "gene_name"
        String sjdb_gtf_tag_exon_parent_gene_type = "gene_type gene_biotype"
        Boolean use_all_cores = false
        Int genome_chr_bin_n_bits = 18
        #@ except: SnakeCase
        Int genome_SA_index_n_bases = 14
        #@ except: SnakeCase
        Int genome_SA_sparse_d = 1
        Int genome_suffix_length_max = -1
        Int sjdb_overhang = 125
        Int ncpu = 8
        Int memory_gb = 50
        Int modify_disk_size_gb = 0
    }

    String star_db_tar_gz = db_name + ".tar.gz"

    Float reference_fasta_size = size(reference_fasta, "GiB")
    Float gtf_size = size(gtf, "GiB")
    Int disk_size_gb = (
        ceil((reference_fasta_size + gtf_size) * 3) + 10 + modify_disk_size_gb
    )

    # Leave 2GB as system overhead
    String memory_limit_bytes = "~{memory_gb - 2}000000000"

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi

        gtf_name=~{basename(gtf, ".gz")}
        gunzip -c ~{gtf} > "$gtf_name" || ln -sf ~{gtf} "$gtf_name"

        ref_fasta=~{basename(reference_fasta, ".gz")}
        gunzip -c ~{reference_fasta} > "$ref_fasta" \
            || ln -sf ~{reference_fasta} "$ref_fasta"

        mkdir ~{db_name};
        STAR --runMode genomeGenerate \
            --genomeDir ~{db_name} \
            --runThreadN "$n_cores" \
            --limitGenomeGenerateRAM ~{memory_limit_bytes} \
            --genomeFastaFiles "$ref_fasta" \
            --sjdbGTFfile "$gtf_name" \
            --sjdbGTFchrPrefix ~{sjdb_gtf_chr_prefix} \
            --sjdbGTFfeatureExon ~{sjdb_gtf_feature_exon} \
            --sjdbGTFtagExonParentTranscript ~{sjdb_gtf_tag_exon_parant_transcript} \
            --sjdbGTFtagExonParentGene ~{sjdb_gtf_tag_exon_parent_gene} \
            --sjdbGTFtagExonParentGeneName ~{sjdb_gtf_tag_exon_parent_gene_name} \
            --sjdbGTFtagExonParentGeneType ~{sjdb_gtf_tag_exon_parent_gene_type} \
            --genomeChrBinNbits ~{genome_chr_bin_n_bits} \
            --genomeSAindexNbases ~{genome_SA_index_n_bases} \
            --genomeSAsparseD ~{genome_SA_sparse_d} \
            --genomeSuffixLengthMax ~{genome_suffix_length_max} \
            --sjdbOverhang ~{sjdb_overhang}

        rm "$gtf_name" "$ref_fasta"

        tar -C ~{db_name} -czf ~{star_db_tar_gz} .
    >>>

    output {
        File star_db = star_db_tar_gz
    }

    runtime {
        cpu: ncpu
        memory: "~{memory_gb} GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/star:2.7.11b--h43eeafb_3"
        maxRetries: 1
    }
}

task alignment {
    meta {
        description: "Runs the STAR aligner on a set of RNA-Seq FASTQ files"
        external_help: "https://github.com/alexdobin/STAR/blob/2.7.11b/doc/STARmanual.pdf"
        outputs: {
            star_log: "Summary mapping statistics after mapping job is complete. The statistics are calculated for each read (Single- or Paired-End) and then summed or averaged over all reads. Note that STAR counts a Paired-End read as one read. Most of the information is collected about the UNIQUE mappers. Each splicing is counted in the numbers of splices, which would correspond to summing the counts in SJ.out.tab. The mismatch/indel error rates are calculated on a per base basis, i.e. as total number of mismatches/indels in all unique mappers divided by the total number of mapped bases.",
            star_bam: "STAR aligned BAM",
            star_junctions: "File contains high confidence collapsed splice junctions in tab-delimited format. Note that STAR defines the junction start/end as intronic bases, while many other software define them as exonic bases. See `meta.external_help` for file specification.",
            star_chimeric_junctions: "Tab delimited file containing chimeric reads and associated metadata. See `meta.external_help` for file specification.",
        }
    }

    parameter_meta {
        star_db_tar_gz: "A gzipped TAR file containing the STAR reference files. The name of the root directory which was archived must match the archive's filename without the `.tar.gz` extension."
        read_one_fastqs_gz: "An array of gzipped FASTQ files containing read one information"
        read_groups: {
            description: "An array of `String`s where each `String` corresponds to one read group.",
            help: "Each read group string should start with the `ID` field followed by any other read group fields, where fields are delimited by a space. See `../data_structures/read_group.wdl` for information about possible fields and utility tasks for constructing, validating, and \"stringifying\" read groups.",
            warning: "The `ID` field for each read group _must_ be contained in the basename of a FASTQ file or pair of FASTQ files if Paired-End. Example: `[\"ID:rg1 PU:flowcell1.lane1 SM:sample1 PL:illumina LB:sample1_lib1\", \"ID:rg2 PU:flowcell1.lane2 SM:sample1 PL:illumina LB:sample1_lib1\"]`. These two read groups could be associated with the following four FASTQs: `[\"sample1.rg1.R1.fastq\", \"sample1.rg2.R1.fastq\"]` and `[\"sample1.rg1.R2.fastq\", \"sample1.rg2.R2.fastq\"]`",
        }
        read_two_fastqs_gz: {
            description: "An array of gzipped FASTQ files containing read two information",
            group: "Common",
        }
        out_sj_filter_intron_max_vs_read_n: {
            description: "maximum gap allowed for junctions supported by 1,2,3,,,N reads. i.e. by default junctions supported by 1 read can have gaps <=50000b, by 2 reads: <=100000b, by 3 reads: <=200000b. by >=4 reads any gap <=alignIntronMax. Does not apply to annotated junctions.",
            group: "Splice Junctions",
        }
        out_sj_filter_overhang_min: {
            description: "minimum overhang length for splice junctions on both sides for: (1) non-canonical motifs, (2) GT/AG and CT/AC motif, (3) GC/AG and CT/GC motif, (4) AT/AC and GT/AT motif. -1 means no output for that motif. Does not apply to annotated junctions.",
            group: "Splice Junctions",
        }
        out_sj_filter_count_unique_min: {
            description: "minimum uniquely mapping read count per junction for: (1) non-canonical motifs, (2) GT/AG and CT/AC motif, (3) GC/AG and CT/GC motif, (4) AT/AC and GT/AT motif. -1 means no output for that motif. Junctions are output if one of outSJfilterCountUniqueMin *OR* outSJfilterCountTotalMin conditions are satisfied. Does not apply to annotated junctions.",
            group: "Splice Junctions",
        }
        out_sj_filter_count_total_min: {
            description: "minimum total (multi-mapping+unique) read count per junction for: (1) non-canonical motifs, (2) GT/AG and CT/AC motif, (3) GC/AG and CT/GC motif, (4) AT/AC and GT/AT motif. -1 means no output for that motif. Junctions are output if one of outSJfilterCountUniqueMin *OR* outSJfilterCountTotalMin conditions are satisfied. Does not apply to annotated junctions.",
            group: "Splice Junctions",
        }
        out_sj_filter_dist_to_other_sj_min: {
            description: "minimum allowed distance to other junctions' donor/acceptor for: (1) non-canonical motifs, (2) GT/AG and CT/AC motif, (3) GC/AG and CT/GC motif, (4) AT/AC and GT/AT motif. Does not apply to annotated junctions.",
            group: "Splice Junctions",
        }
        align_sj_stitch_mismatch_n_max: {
            description: "maximum number of mismatches for stitching of the splice junctions (-1: no limit) for: (1) non-canonical motifs, (2) GT/AG and CT/AC motif, (3) GC/AG and CT/GC motif, (4) AT/AC and GT/AT motif",
            group: "Splice Junctions",
        }
        clip_3p_adapter_seq: {
            description: "adapter sequences to clip from 3p of each mate. `left` applies to read one and `right` applies to read two.",
            choices: {
                None: "No 3p adapter trimming will be performed",
                sequence: "A nucleotide sequence string of any length, matching the regex `/[ATCG]+/`",
                polyA: "polyA sequence with the length equal to read length",
            },
            group: "Common",
        }
        clip_3p_adapter_mmp: "max proportion of mismatches for 3p adapter clipping for each mate. `left` applies to read one and `right` applies to read two."
        align_ends_protrude: {
            description: "allow protrusion of alignment ends, i.e. start (end) of the +strand mate downstream of the start (end) of the -strand mate. `left`: maximum number of protrusion bases allowed. `right`: see `choices` below.",
            choices: {
                ConcordantPair: "report alignments with non-zero protrusion as concordant pairs",
                DiscordantPair: "report alignments with non-zero protrusion as discordant pairs",
            },
        }
        clip_3p_n_bases: "number of bases to clip from 3p of each mate. `left` applies to read one and `right` applies to read two."
        clip_3p_after_adapter_n_bases: "number of bases to clip from 3p of each mate after the adapter clipping. `left` applies to read one and `right` applies to read two."
        clip_5p_n_bases: "number of bases to clip from 5p of each mate. `left` applies to read one and `right` applies to read two."
        prefix: {
            description: "Prefix for the BAM and other STAR files. The extensions `.Aligned.out.bam`, `.Log.final.out`, `.SJ.out.tab`, and `.Chimeric.out.junction` will be added.",
            help: "See `../README.md` for more information on the default prefix evaluation.",
            group: "Common",
        }
        read_name_separator: "character(s) separating the part of the read names that will be trimmed in output (read name after space is always trimmed)"
        clip_adapter_type: {
            description: "adapter clipping type",
            choices: {
                Hamming: "adapter clipping based on Hamming distance, with the number of mismatches controlled by --clip5pAdapterMMp",
                CellRanger4: "5p and 3p adapter clipping similar to CellRanger4. Utilizes Opal package by Martin Šošić: https://github.com/Martinsos/opal",
                None: "no adapter clipping, all other clip* parameters are disregarded",
            },
        }
        out_sam_strand_field: {
            description: "Cufflinks-like strand field flag",
            choices: {
                None: "not used",
                intronMotif: "strand derived from the intron motif. This option changes the output alignments: reads with inconsistent and/or non-canonical introns are filtered out.",
            },
            group: "Common",
        }
        out_sam_attributes: {
            description: "a string of desired SAM attributes, in the order desired for the output SAM. Tags can be listed in any combination/order. **[STAR defaults]**: `NH HI AS nM`. **[WDL default]**: `NH HI AS nM NM MD XS`.",
            choices: {
                NH: "number of loci the reads maps to: =1 for unique mappers, >1 for multimappers. Standard SAM tag.",
                HI: "multiple alignment index, starts with --outSAMattrIHstart (=1 by default). Standard SAM tag.",
                AS: "local alignment score, +1/-1 for matches/mismateches, score* penalties for indels and gaps. For PE reads, total score for two mates. Standard SAM tag.",
                nM: "number of mismatches. For PE reads, sum over two mates.",
                NM: "edit distance to the reference (number of mismatched + inserted + deleted bases) for each mate. Standard SAM tag.",
                MD: "string encoding mismatched and deleted reference bases (see standard SAM specifications). Standard SAM tag.",
                jM: "intron motifs for all junctions (i.e. N in CIGAR): 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5: AT/AC, 6: GT/AT. If splice junctions database is used, and a junction is annotated, 20 is added to its motif value.",
                jI: "start and end of introns for all junctions (1-based).",
                XS: "alignment strand according to --outSAMstrandField.",
                MC: "mate's CIGAR string. Standard SAM tag.",
                ch: "marks all segments of all chimeric alignments for --chimOutType WithinBAM output.",
                cN: "number of bases clipped from the read ends: 5' and 3'",
            },
            group: "Common",
        }
        out_sam_unmapped: {
            description: "output of unmapped reads in the SAM format.",
            choices: {
                None: "no output **[STAR default]**",
                Within: "output unmapped reads within the main SAM file (i.e. Aligned.out.sam) **[WDL default]**",
            },
        }
        out_sam_order: {
            description: "type of sorting for the SAM output",
            choices: {
                Paired: "one mate after the other for all paired alignments",
                PairedKeepInputOrder: "one mate after the other for all paired alignments, the order is kept the same as in the input FASTQ files",
            },
        }
        out_sam_read_id: {
            description: "read ID record type",
            choices: {
                Standard: "first word (until space) from the FASTx read ID line, removing /1,/2 from the end",
                Number: "read number (index) in the FASTx file",
            },
        }
        out_sam_tlen: {
            description: "calculation method for the TLEN field in the SAM/BAM files",
            choices: {
                left_plus: "leftmost base of the (+)strand mate to rightmost base of the (-)mate. (+)sign for the (+)strand mate",
                left_any: "leftmost base of any mate to rightmost base of any mate. (+)sign for the mate with the leftmost base. This is different from `left_plus` for overlapping mates with protruding ends",
            },
        }
        out_filter_type: {
            description: "type of filtering",
            choices: {
                Normal: "standard filtering using only current alignment",
                BySJout: "keep only those reads that contain junctions that passed filtering into SJ.out.tab",
            },
            group: "Common",
        }
        out_filter_intron_motifs: {
            description: "filter alignment using their motifs",
            choices: {
                None: "no filtering",
                RemoveNoncanonical: "filter out alignments that contain non-canonical junctions",
                RemoveNoncanonicalUnannotated: "filter out alignments that contain non-canonical unannotated junctions when using annotated splice junctions database. The annotated non-canonical junctions will be kept.",
            },
            group: "Common",
        }
        out_filter_intron_strands: {
            description: "filter alignments",
            choices: {
                None: "no filtering",
                RemoveInconsistentStrands: "remove alignments that have junctions with inconsistent strands",
            },
            group: "Common",
        }
        out_sj_filter_reads: {
            description: "which reads to consider for collapsed splice junctions output",
            choices: {
                All: "all reads, unique- and multi-mappers",
                Unique: "uniquely mapping reads only",
            },
            group: "Common",
        }
        align_ends_type: {
            description: "type of read ends alignment",
            choices: {
                Local: "standard local alignment with soft-clipping allowed",
                EndToEnd: "force end-to-end read alignment, do not soft-clip",
                Extend5pOfRead1: "fully extend only the 5p of the read1, all other ends: local alignment",
                Extend5pOfReads12: "fully extend only the 5p of the both read1 and read2, all other ends: local alignment",
            },
        }
        align_soft_clip_at_reference_ends: {
            description: "allow the soft-clipping of the alignments past the end of the chromosomes",
            choices: {
                Yes: "allow",
                No: "prohibit, useful for compatibility with Cufflinks",
            },
            group: "Common",
        }
        align_insertion_flush: {
            description: "how to flush ambiguous insertion positions",
            choices: {
                None: "insertions are not flushed",
                Right: "insertions are flushed to the right",
            },
            group: "Common",
        }
        chim_out_type: {
            description: "type of chimeric output",
            choices: {
                Junctions: "Chimeric.out.junction",
                WithinBAM_HardClip: "output into main aligned BAM files (Aligned.*.bam). Hard-clipping in the CIGAR for supplemental chimeric alignments.",
                WithinBAM_SoftClip: "output into main aligned BAM files (Aligned.*.bam). Soft-clipping in the CIGAR for supplemental chimeric alignments.",
            },
            tool_default: "Junctions",
            group: "Common",
        }
        chim_filter: {
            description: "different filters for chimeric alignments",
            choices: {
                None: "no filtering",
                banGenomicN: "Ns are not allowed in the genome sequence around the chimeric junction",
            },
        }
        chim_out_junction_format: {
            description: "formatting type for the Chimeric.out.junction file",
            choices: {
                plain: "no comment lines/headers",
                comments: "comment lines at the end of the file: command line and Nreads: total, unique/multi-mapping",
            },
            group: "Common",
        }
        twopass_mode: {
            description: "2-pass mapping mode",
            choices: {
                None: "1-pass mapping **[STAR default]**",
                Basic: "basic 2-pass mapping, with all 1st pass junctions inserted into the genome indices on the fly **[WDL default]**",
            },
            group: "Common",
        }
        use_all_cores: {
            description: "Use all cores? Recommended for cloud environments.",
            group: "Resources",
        }
        out_filter_mismatch_n_over_l_max: "alignment will be output only if its ratio of mismatches to *mapped* length is less than or equal to this value"
        out_filter_mismatch_n_over_read_l_max: "alignment will be output only if its ratio of mismatches to *read* length is less than or equal to this value"
        out_filter_score_min_over_l_read: "same as outFilterScoreMin, but normalized to read length (sum of mates' lengths for Paired-End reads)"
        out_filter_match_n_min_over_l_read: "same as outFilterMatchNmin, but normalized to the read length (sum of mates' lengths for Paired-End reads)"
        score_genomic_length_log2_scale: "extra score logarithmically scaled with genomic length of the alignment: scoreGenomicLengthLog2scale*log2(genomicLength)"
        seed_search_start_l_max_over_l_read: "seedSearchStartLmax normalized to read length (sum of mates' lengths for Paired-End reads)"
        align_spliced_mate_map_l_min_over_l_mate: "alignSplicedMateMapLmin normalized to mate length"
        pe_overlap_mmp: "maximum proportion of mismatched bases in the overlap area"
        run_rng_seed: {
            description: "random number generator seed",
            group: "Common",
        }
        sjdb_score: {
            description: "extra alignment score for alignments that cross database junctions",
            group: "Common",
        }
        read_map_number: {
            description: "number of reads to map from the beginning of the file. -1 to map all reads",
            group: "Common",
        }
        read_quality_score_base: "number to be subtracted from the ASCII code to get Phred quality score"
        limit_out_sj_one_read: "max number of junctions for one read (including all multi-mappers)"
        limit_out_sj_collapsed: "max number of collapsed junctions"
        limit_sjdb_insert_n_sj: "maximum number of junction to be inserted to the genome on the fly at the mapping stage, including those from annotations and those detected in the 1st step of the 2-pass run"
        out_QS_conversion_add: "add this number to the quality score (e.g. to convert from Illumina to Sanger, use -31)"
        out_sam_attr_IH_start: "start value for the IH attribute. 0 may be required by some downstream software, such as Cufflinks or StringTie."
        out_sam_mapq_unique: "`0-255`: the MAPQ value for unique mappers. Please note the STAR default (255) produces errors downstream, as a MAPQ value of 255 is reserved to indicate a missing value. The default of this task is 254, which is the highest _valid_ MAPQ value, and possibly what the author of STAR intended. **[STAR default]**: `255`. **[WDL default]**: `254`."
        out_sam_flag_OR: "`0-65535`: sam FLAG will be bitwise OR'd with this value, i.e. FLAG=FLAG | outSAMflagOR. This is applied after all flags have been set by STAR, and after outSAMflagAND. Can be used to set specific bits that are not set otherwise."
        out_sam_flag_AND: "`0-65535`: sam FLAG will be bitwise AND'd with this value, i.e. FLAG=FLAG & outSAMflagOR. This is applied after all flags have been set by STAR, but before outSAMflagOR. Can be used to unset specific bits that are not set otherwise."
        out_filter_multimap_score_range: "the score range below the maximum score for multimapping alignments"
        out_filter_multimap_n_max: {
            description: "maximum number of loci the read is allowed to map to. Alignments (all of them) will be output only if the read maps to no more loci than this value. Otherwise no alignments will be output, and the read will be counted as 'mapped to too many loci' in the Log.final.out. **[STAR default]**: `10`. **[WDL default]**: `20`.",
            group: "Common",
        }
        out_filter_mismatch_n_max: {
            description: "alignment will be output only if it has no more mismatches than this value",
            group: "Common",
        }
        out_filter_score_min: {
            description: "alignment will be output only if its score is higher than or equal to this value",
            group: "Common",
        }
        out_filter_match_n_min: {
            description: "alignment will be output only if the number of matched bases is higher than or equal to this value",
            group: "Common",
        }
        score_gap: "splice junction penalty (independent on intron motif)"
        score_gap_noncanon: "non-canonical junction penalty (in addition to scoreGap)"
        score_gap_GCAG: "GC/AG and CT/GC junction penalty (in addition to scoreGap)"
        score_gap_ATAC: "AT/AC and GT/AT junction penalty (in addition to scoreGap)"
        score_del_open: "deletion open penalty"
        score_del_base: "deletion extension penalty per base (in addition to scoreDelOpen)"
        score_ins_open: "insertion open penalty"
        score_ins_base: "insertion extension penalty per base (in addition to scoreInsOpen)"
        score_stitch_sj_shift: "maximum score reduction while searching for SJ boundaries in the stitching step"
        seed_search_start_l_max: "defines the search start point through the read - the read is split into pieces no longer than this value"
        seed_search_l_max: "defines the maximum length of the seeds, if =0 seed length is not limited"
        seed_multimap_n_max: "only pieces that map fewer than this value are utilized in the stitching procedure"
        seed_per_read_n_max: "max number of seeds per read"
        seed_per_window_n_max: "max number of seeds per window"
        seed_none_loci_per_window: "max number of one seed loci per window"
        seed_split_min: "min length of the seed sequences split by Ns or mate gap"
        seed_map_min: "min length of seeds to be mapped"
        align_intron_min: {
            description: "minimum intron size: genomic gap is considered intron if its length>=alignIntronMin, otherwise it is considered Deletion",
            group: "Common",
        }
        align_intron_max: {
            description: "maximum intron size, if 0, max intron size will be determined by (2^winBinNbits)*winAnchorDistNbins. **[STAR default]**: `0`. **[WDL default]**: `500000`.",
            group: "Common",
        }
        align_mates_gap_max: {
            description: "maximum gap between two mates, if 0, max intron gap will be determined by (2^winBinNbits)*winAnchorDistNbins. **[STAR default]**: `0`. **[WDL default]**: `1000000`",
            group: "Common",
        }
        align_sj_overhang_min: {
            description: "minimum overhang (i.e. block size) for spliced alignments",
            group: "Common",
        }
        align_sjdb_overhang_min: {
            description: "minimum overhang (i.e. block size) for annotated (sjdb) spliced alignments. **[STAR default]**: `3`. **[WDL default]**: `1`.",
            group: "Common",
        }
        align_spliced_mate_map_l_min: "minimum mapped length for a read mate that is spliced"
        align_windows_per_read_n_max: "max number of windows per read"
        align_transcripts_per_window_n_max: "max number of transcripts per window"
        align_transcripts_per_read_n_max: "max number of different alignments per read to consider"
        pe_overlap_n_bases_min: "minimum number of overlap bases to trigger mates merging and realignment. Specify >0 value to switch on the 'merging of overlapping mates' algorithm."
        win_anchor_multimap_n_max: "max number of loci anchors are allowed to map to"
        win_bin_n_bits: "=log2(winBin), where winBin is the size of the bin for the windows/clustering, each window will occupy an integer number of bins"
        win_anchor_dist_n_bins: "max number of bins between two anchors that allows aggregation of anchors into one window"
        win_flank_n_bins: "=log2(winFlank), where winFlank is the size of the left and right flanking regions for each window"
        chim_segment_min: {
            description: "minimum length of chimeric segment length, if ==0, no chimeric output",
            tool_default: 0,
            group: "Common",
        }
        chim_score_min: {
            description: "minimum total (summed) score of the chimeric segments",
            group: "Common",
        }
        chim_score_drop_max: {
            description: "max drop (difference) of chimeric score (the sum of scores of all chimeric segments) from the read length",
            group: "Common",
        }
        chim_score_separation: "minimum difference (separation) between the best chimeric score and the next one"
        chim_score_junction_nonGTAG: "penalty for a non-GT/AG chimeric junction"
        chim_junction_overhang_min: {
            description: "minimum overhang for a chimeric junction",
            group: "Common",
        }
        chim_segment_read_gap_max: {
            description: "maximum gap in the read sequence between chimeric segments",
            group: "Common",
        }
        chim_main_segment_multi_n_max: {
            description: "maximum number of multi-alignments for the main chimeric segment. =1 will prohibit multimapping main segments.",
            group: "Common",
        }
        chim_multimap_n_max: {
            description: "maximum number of chimeric multi-alignments. `0`: use the old scheme for chimeric detection which only considered unique alignments",
            group: "Common",
        }
        chim_multimap_score_range: "the score range for multi-mapping chimeras below the best chimeric score. Only works with --chimMultimapNmax > 1."
        chim_nonchim_score_drop_min: "to trigger chimeric detection, the drop in the best non-chimeric alignment score with respect to the read length has to be greater than this value"
        twopass1_reads_n: {
            description: "number of reads to process for the 1st step. Use default (`-1`) to map all reads in the first step",
            group: "Common",
        }
        ncpu: {
            description: "Number of cores to allocate for task",
            group: "Resources",
        }
        modify_disk_size_gb: {
            description: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB.",
            group: "Resources",
        }
    }

    input {
        File star_db_tar_gz
        Array[File]+ read_one_fastqs_gz
        Array[String]+ read_groups
        Array[File] read_two_fastqs_gz = []
        Array[Int] out_sj_filter_intron_max_vs_read_n = [50000, 100000, 200000]
        SpliceJunctionMotifs out_sj_filter_overhang_min = SpliceJunctionMotifs {
            noncanonical_motifs: 30,
            GT_AG_and_CT_AC_motif: 12,
            GC_AG_and_CT_GC_motif: 12,
            AT_AC_and_GT_AT_motif: 12,
        }
        SpliceJunctionMotifs out_sj_filter_count_unique_min = SpliceJunctionMotifs {
            noncanonical_motifs: 3,
            GT_AG_and_CT_AC_motif: 1,
            GC_AG_and_CT_GC_motif: 1,
            AT_AC_and_GT_AT_motif: 1,
        }
        SpliceJunctionMotifs out_sj_filter_count_total_min = SpliceJunctionMotifs {
            noncanonical_motifs: 3,
            GT_AG_and_CT_AC_motif: 1,
            GC_AG_and_CT_GC_motif: 1,
            AT_AC_and_GT_AT_motif: 1,
        }
        SpliceJunctionMotifs out_sj_filter_dist_to_other_sj_min = SpliceJunctionMotifs {
            noncanonical_motifs: 10,
            GT_AG_and_CT_AC_motif: 0,
            GC_AG_and_CT_GC_motif: 5,
            AT_AC_and_GT_AT_motif: 10,
        }
        SpliceJunctionMotifs align_sj_stitch_mismatch_n_max = SpliceJunctionMotifs {
            noncanonical_motifs: 0,
            GT_AG_and_CT_AC_motif: -1,
            GC_AG_and_CT_GC_motif: 0,
            AT_AC_and_GT_AT_motif: 0,
        }
        Pair[String, String] clip_3p_adapter_seq = ("None", "None")
        Pair[Float, Float] clip_3p_adapter_mmp = (0.1, 0.1)
        Pair[Int, String] align_ends_protrude = (0, "ConcordantPair")
        Pair[Int, Int] clip_3p_n_bases = (0, 0)
        Pair[Int, Int] clip_3p_after_adapter_n_bases = (0, 0)
        Pair[Int, Int] clip_5p_n_bases = (0, 0)
        String prefix = sub(
            basename(read_one_fastqs_gz[0]),
            "(([_.][rR](?:ead)?[12])((?:[_.-][^_.-]*?)*?))?\\.(fastq|fq)(\\.gz)?$",
            ""  # Once replacing with capturing groups is supported, replace with group 3
        )
        String read_name_separator = "/"
        String clip_adapter_type = "Hamming"
        String out_sam_strand_field = "intronMotif"
        String out_sam_attributes = "NH HI AS nM NM MD XS"
        String out_sam_unmapped = "Within"
        String out_sam_order = "Paired"
        String out_sam_read_id = "Standard"
        String out_sam_tlen = "left_plus"
        String out_filter_type = "Normal"
        String out_filter_intron_motifs = "None"
        String out_filter_intron_strands = "RemoveInconsistentStrands"
        String out_sj_filter_reads = "All"
        String align_ends_type = "Local"
        String align_soft_clip_at_reference_ends = "Yes"
        String align_insertion_flush = "None"
        String chim_out_type = "WithinBAM HardClip"
        String chim_filter = "banGenomicN"
        String chim_out_junction_format = "plain"
        String twopass_mode = "Basic"
        Boolean use_all_cores = false
        Float out_filter_mismatch_n_over_l_max = 0.3
        Float out_filter_mismatch_n_over_read_l_max = 1.0
        Float out_filter_score_min_over_l_read = 0.66
        Float out_filter_match_n_min_over_l_read = 0.66
        Float score_genomic_length_log2_scale = -0.25
        Float seed_search_start_l_max_over_l_read = 1.0
        Float align_spliced_mate_map_l_min_over_l_mate = 0.66
        Float pe_overlap_mmp = 0.01
        Int run_rng_seed = 777
        Int sjdb_score = 2
        Int read_map_number = -1
        Int read_quality_score_base = 33
        Int limit_out_sj_one_read = 1000
        Int limit_out_sj_collapsed = 1000000
        Int limit_sjdb_insert_n_sj = 1000000
        #@ except: SnakeCase
        Int out_QS_conversion_add = 0
        #@ except: SnakeCase
        Int out_sam_attr_IH_start = 1
        Int out_sam_mapq_unique = 254
        #@ except: SnakeCase
        Int out_sam_flag_OR = 0
        #@ except: SnakeCase
        Int out_sam_flag_AND = 65535
        Int out_filter_multimap_score_range = 1
        Int out_filter_multimap_n_max = 10
        Int out_filter_mismatch_n_max = 10
        Int out_filter_score_min = 0
        Int out_filter_match_n_min = 0
        Int score_gap = 0
        Int score_gap_noncanon = -8
        #@ except: SnakeCase
        Int score_gap_GCAG = -4
        #@ except: SnakeCase
        Int score_gap_ATAC = -8
        Int score_del_open = -2
        Int score_del_base = -2
        Int score_ins_open = -2
        Int score_ins_base = -2
        Int score_stitch_sj_shift = 1
        Int seed_search_start_l_max = 50
        Int seed_search_l_max = 0
        Int seed_multimap_n_max = 10000
        Int seed_per_read_n_max = 1000
        Int seed_per_window_n_max = 50
        Int seed_none_loci_per_window = 10
        Int seed_split_min = 12
        Int seed_map_min = 5
        Int align_intron_min = 21
        Int align_intron_max = 500000
        Int align_mates_gap_max = 1000000
        Int align_sj_overhang_min = 5
        Int align_sjdb_overhang_min = 1
        Int align_spliced_mate_map_l_min = 0
        Int align_windows_per_read_n_max = 10000
        Int align_transcripts_per_window_n_max = 100
        Int align_transcripts_per_read_n_max = 10000
        Int pe_overlap_n_bases_min = 0
        Int win_anchor_multimap_n_max = 50
        Int win_bin_n_bits = 16
        Int win_anchor_dist_n_bins = 9
        Int win_flank_n_bins = 4
        Int chim_segment_min = 10
        Int chim_score_min = 0
        Int chim_score_drop_max = 20
        Int chim_score_separation = 10
        #@ except: SnakeCase
        Int chim_score_junction_nonGTAG = -1
        Int chim_junction_overhang_min = 20
        Int chim_segment_read_gap_max = 0
        Int chim_main_segment_multi_n_max = 10
        Int chim_multimap_n_max = 0
        Int chim_multimap_score_range = 1
        Int chim_nonchim_score_drop_min = 20
        Int twopass1_reads_n = -1
        Int ncpu = 8
        Int modify_disk_size_gb = 0
    }

    Float read_one_fastqs_size = size(read_one_fastqs_gz, "GiB")
    Float read_two_fastqs_size = size(read_two_fastqs_gz, "GiB")
    Float star_db_tar_gz_size = size(star_db_tar_gz, "GiB")
    Int disk_size_gb = (
        (
            ceil(read_one_fastqs_size + read_two_fastqs_size + star_db_tar_gz_size) * 3
        ) + 10 + modify_disk_size_gb
    )

    #@ except: LineWidth
    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi

        mkdir star_db
        tar -xzf ~{star_db_tar_gz} -C star_db/ --no-same-owner

        STAR --readFilesIn ~{sep(",", squote(read_one_fastqs_gz))} ~{sep(",", squote(read_two_fastqs_gz))} \
            --readFilesCommand "gunzip -c" \
            --genomeDir star_db \
            --runThreadN "$n_cores" \
            --outSAMtype BAM Unsorted \
            --outMultimapperOrder Random \
            --outFileNamePrefix ~{prefix + "."} \
            --twopassMode ~{twopass_mode} \
            --outSAMattrRGline ~{sep(" , ", squote(read_groups))} \
            --outSJfilterIntronMaxVsReadN ~{
                sep(" ", quote(out_sj_filter_intron_max_vs_read_n))
            } \
            --outSJfilterOverhangMin ~{sep(" ", quote([
                out_sj_filter_overhang_min.noncanonical_motifs,
                out_sj_filter_overhang_min.GT_AG_and_CT_AC_motif,
                out_sj_filter_overhang_min.GC_AG_and_CT_GC_motif,
                out_sj_filter_overhang_min.AT_AC_and_GT_AT_motif,
            ]))} \
            --outSJfilterCountUniqueMin ~{sep(" ", quote([
                out_sj_filter_count_unique_min.noncanonical_motifs,
                out_sj_filter_count_unique_min.GT_AG_and_CT_AC_motif,
                out_sj_filter_count_unique_min.GC_AG_and_CT_GC_motif,
                out_sj_filter_count_unique_min.AT_AC_and_GT_AT_motif,
            ]))} \
            --outSJfilterCountTotalMin ~{sep(" ", quote([
                out_sj_filter_count_total_min.noncanonical_motifs,
                out_sj_filter_count_total_min.GT_AG_and_CT_AC_motif,
                out_sj_filter_count_total_min.GC_AG_and_CT_GC_motif,
                out_sj_filter_count_total_min.AT_AC_and_GT_AT_motif,
            ]))} \
            --outSJfilterDistToOtherSJmin ~{sep(" ", quote([
                out_sj_filter_dist_to_other_sj_min.noncanonical_motifs,
                out_sj_filter_dist_to_other_sj_min.GT_AG_and_CT_AC_motif,
                out_sj_filter_dist_to_other_sj_min.GC_AG_and_CT_GC_motif,
                out_sj_filter_dist_to_other_sj_min.AT_AC_and_GT_AT_motif,
            ]))} \
            --alignSJstitchMismatchNmax ~{sep(" ", quote([
                align_sj_stitch_mismatch_n_max.noncanonical_motifs,
                align_sj_stitch_mismatch_n_max.GT_AG_and_CT_AC_motif,
                align_sj_stitch_mismatch_n_max.GC_AG_and_CT_GC_motif,
                align_sj_stitch_mismatch_n_max.AT_AC_and_GT_AT_motif,
            ]))} \
            --clip3pAdapterSeq ~{clip_3p_adapter_seq.left} ~{(
                if (length(read_two_fastqs_gz) != 0)
                then clip_3p_adapter_seq.right
                else ""
            )} \
            --clip3pAdapterMMp ~{clip_3p_adapter_mmp.left} ~{(
                if (length(read_two_fastqs_gz) != 0)
                then clip_3p_adapter_mmp.right
                else None
            )} \
            --alignEndsProtrude ~{align_ends_protrude.left} ~{(
                if (length(read_two_fastqs_gz) != 0)
                then align_ends_protrude.right
                else None
            )} \
            --clip3pNbases ~{clip_3p_n_bases.left} ~{(
                if (length(read_two_fastqs_gz) != 0)
                then clip_3p_n_bases.right
                else None
            )} \
            --clip3pAfterAdapterNbases ~{clip_3p_after_adapter_n_bases.left} ~{(
                if (length(read_two_fastqs_gz) != 0)
                then clip_3p_after_adapter_n_bases.right
                else None
            )} \
            --clip5pNbases ~{clip_5p_n_bases.left} ~{(
                if (length(read_two_fastqs_gz) != 0)
                then clip_5p_n_bases.right
                else None
            )} \
            --readNameSeparator ~{read_name_separator} \
            --clipAdapterType ~{clip_adapter_type} \
            --outSAMstrandField ~{out_sam_strand_field} \
            --outSAMattributes ~{out_sam_attributes} \
            --outSAMunmapped ~{out_sam_unmapped} \
            --outSAMorder ~{out_sam_order} \
            --outSAMreadID ~{out_sam_read_id} \
            --outSAMtlen ~{(
                if (out_sam_tlen == "left_plus")
                then "1"
                else (
                    if (out_sam_tlen == "left_any") then "2" else "error"
                )
            )} \
            --outFilterType ~{out_filter_type} \
            --outFilterIntronMotifs ~{out_filter_intron_motifs} \
            --outFilterIntronStrands ~{out_filter_intron_strands} \
            --outSJfilterReads ~{out_sj_filter_reads} \
            --alignEndsType ~{align_ends_type} \
            --alignSoftClipAtReferenceEnds ~{align_soft_clip_at_reference_ends} \
            --alignInsertionFlush ~{align_insertion_flush} \
            --chimOutType ~{chim_out_type} \
            --chimFilter ~{chim_filter} \
            --chimOutJunctionFormat ~{chim_out_junction_format} \
            --outFilterMismatchNoverLmax ~{out_filter_mismatch_n_over_l_max} \
            --outFilterMismatchNoverReadLmax ~{out_filter_mismatch_n_over_read_l_max} \
            --outFilterScoreMinOverLread ~{out_filter_score_min_over_l_read} \
            --outFilterMatchNminOverLread ~{out_filter_match_n_min_over_l_read} \
            --scoreGenomicLengthLog2scale ~{score_genomic_length_log2_scale} \
            --seedSearchStartLmaxOverLread ~{seed_search_start_l_max_over_l_read} \
            --alignSplicedMateMapLminOverLmate ~{align_spliced_mate_map_l_min_over_l_mate} \
            --peOverlapMMp ~{pe_overlap_mmp} \
            --runRNGseed ~{run_rng_seed} \
            --sjdbScore ~{sjdb_score} \
            --readMapNumber ~{read_map_number} \
            --readQualityScoreBase ~{read_quality_score_base} \
            --limitOutSJoneRead ~{limit_out_sj_one_read} \
            --limitOutSJcollapsed ~{limit_out_sj_collapsed} \
            --limitSjdbInsertNsj ~{limit_sjdb_insert_n_sj} \
            --outQSconversionAdd ~{out_QS_conversion_add} \
            --outSAMattrIHstart ~{out_sam_attr_IH_start} \
            --outSAMmapqUnique ~{out_sam_mapq_unique} \
            --outSAMflagOR ~{out_sam_flag_OR} \
            --outSAMflagAND ~{out_sam_flag_AND} \
            --outFilterMultimapScoreRange ~{out_filter_multimap_score_range} \
            --outFilterMultimapNmax ~{out_filter_multimap_n_max} \
            --outFilterMismatchNmax ~{out_filter_mismatch_n_max} \
            --outFilterScoreMin ~{out_filter_score_min} \
            --outFilterMatchNmin ~{out_filter_match_n_min} \
            --scoreGap ~{score_gap} \
            --scoreGapNoncan ~{score_gap_noncanon} \
            --scoreGapGCAG ~{score_gap_GCAG} \
            --scoreGapATAC ~{score_gap_ATAC} \
            --scoreDelOpen ~{score_del_open} \
            --scoreDelBase ~{score_del_base} \
            --scoreInsOpen ~{score_ins_open} \
            --scoreInsBase ~{score_ins_base} \
            --scoreStitchSJshift ~{score_stitch_sj_shift} \
            --seedSearchStartLmax ~{seed_search_start_l_max} \
            --seedSearchLmax ~{seed_search_l_max} \
            --seedMultimapNmax ~{seed_multimap_n_max} \
            --seedPerReadNmax ~{seed_per_read_n_max} \
            --seedPerWindowNmax ~{seed_per_window_n_max} \
            --seedNoneLociPerWindow ~{seed_none_loci_per_window} \
            --seedSplitMin ~{seed_split_min} \
            --seedMapMin ~{seed_map_min} \
            --alignIntronMin ~{align_intron_min} \
            --alignIntronMax ~{align_intron_max} \
            --alignMatesGapMax ~{align_mates_gap_max} \
            --alignSJoverhangMin ~{align_sj_overhang_min} \
            --alignSJDBoverhangMin ~{align_sjdb_overhang_min} \
            --alignSplicedMateMapLmin ~{align_spliced_mate_map_l_min} \
            --alignWindowsPerReadNmax ~{align_windows_per_read_n_max} \
            --alignTranscriptsPerWindowNmax ~{align_transcripts_per_window_n_max} \
            --alignTranscriptsPerReadNmax ~{align_transcripts_per_read_n_max} \
            --peOverlapNbasesMin ~{pe_overlap_n_bases_min} \
            --winAnchorMultimapNmax ~{win_anchor_multimap_n_max} \
            --winBinNbits ~{win_bin_n_bits} \
            --winAnchorDistNbins ~{win_anchor_dist_n_bins} \
            --winFlankNbins ~{win_flank_n_bins} \
            --chimSegmentMin ~{chim_segment_min} \
            --chimScoreMin ~{chim_score_min} \
            --chimScoreDropMax ~{chim_score_drop_max} \
            --chimScoreSeparation ~{chim_score_separation} \
            --chimScoreJunctionNonGTAG ~{chim_score_junction_nonGTAG} \
            --chimJunctionOverhangMin ~{chim_junction_overhang_min} \
            --chimSegmentReadGapMax ~{chim_segment_read_gap_max} \
            --chimMainSegmentMultNmax ~{chim_main_segment_multi_n_max} \
            --chimMultimapNmax ~{chim_multimap_n_max} \
            --chimMultimapScoreRange ~{chim_multimap_score_range} \
            --chimNonchimScoreDropMin ~{chim_nonchim_score_drop_min} \
            --twopass1readsN ~{twopass1_reads_n}
    >>>

    output {
        File star_log = prefix + ".Log.final.out"
        File star_bam = prefix + ".Aligned.out.bam"
        File star_junctions = prefix + ".SJ.out.tab"
        File? star_chimeric_junctions = prefix + ".Chimeric.out.junction"
    }

    runtime {
        cpu: ncpu
        memory: "50 GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/star:2.7.11b--h43eeafb_3"
        maxRetries: 1
    }
}

# There are multiple Splice Junction Motif arguments for STAR
# that are all formatted the same. Use this struct for consistency.
#@ except: SnakeCase
struct SpliceJunctionMotifs {
    Int noncanonical_motifs
    Int GT_AG_and_CT_AC_motif
    Int GC_AG_and_CT_GC_motif
    Int AT_AC_and_GT_AT_motif
}
