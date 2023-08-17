## # STAR
##
## This WDL file wraps the [STAR aligner](https://github.com/alexdobin/STAR).
## STAR is an RNA-Seq aligner.

version 1.1

task build_star_db {
    meta {
        description: "This WDL task runs STAR's build command to generate a STAR format reference for alignment." 
    }

    parameter_meta {
        reference_fasta: "The FASTA format reference file for the genome"
        gtf: "GTF format feature file"
        db_name: "Name for output in compressed, archived format. The suffix `.tar.gz` will be added."
        sjdbGTFchrPrefix: "prefix for chromosome names in a GTF file (e.g. 'chr' for using ENSMEBL annotations with UCSC genomes)"
        sjdbGTFfeatureExon: "feature type in GTF file to be used as exons for building transcripts"
        sjdbGTFtagExonParentTranscript: "GTF attribute name for parent transcript ID"
        sjdbGTFtagExonParentGene: "GTF attribute name for parent gene ID"
        sjdbGTFtagExonParentGeneName: "GTF attrbute name for parent gene name"
        sjdbGTFtagExonParentGeneType: "GTF attrbute name for parent gene type"
        use_all_cores: "Use all cores? Recommended for cloud environments. Not recommended for cluster environments."
        genomeChrBinNbits: "=log2(chrBin), where chrBin is the size of the bins for genome storage: each chromosome will occupy an integer number of bins. For a genome with large number of contigs, it is recommended to scale this parameter as min(18, log2[max(GenomeLength/NumberOfReferences,ReadLength)])."
        genomeSAindexNbases: "length (bases) of the SA pre-indexing string. Typically between 10 and 15. Longer strings will use much more memory, but allow faster searches. For small genomes, the parameter `--genomeSAindexNbases` must be scaled down to `min(14, log2(GenomeLength)/2 - 1)`."
        genomeSAsparseD: "suffix array sparsity, i.e. distance between indices: use bigger numbers to decrease needed RAM at the cost of mapping speed reduction."
        genomeSuffixLengthMax: "maximum length of the suffixes, has to be longer than read length. -1 = infinite."
        sjdbOverhang: "length of the donor/acceptor sequence on each side of the junctions, ideally = (mate_length - 1). **[STAR default]**: `100`. **[WDL default]**: `125`."
        ncpu: "Number of cores to allocate for task"
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
        max_retries: "Number of times to retry in case of failure"
    }

    input {
        File reference_fasta
        File gtf
        String db_name = "star_db"
        String sjdbGTFchrPrefix = "-"
        String sjdbGTFfeatureExon = "exon"
        String sjdbGTFtagExonParentTranscript = "transcript_id"
        String sjdbGTFtagExonParentGene = "gene_id"
        String sjdbGTFtagExonParentGeneName = "gene_name"
        String sjdbGTFtagExonParentGeneType = "gene_type gene_biotype"
        Boolean use_all_cores = false
        Int genomeChrBinNbits = 18
        Int genomeSAindexNbases = 14
        Int genomeSAsparseD = 1
        Int genomeSuffixLengthMax = -1
        Int sjdbOverhang = 125
        Int ncpu = 1
        Int memory_gb = 50
        Int modify_disk_size_gb = 0
        Int max_retries = 1
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
            --sjdbGTFchrPrefix ~{sjdbGTFchrPrefix} \
            --sjdbGTFfeatureExon ~{sjdbGTFfeatureExon} \
            --sjdbGTFtagExonParentTranscript ~{sjdbGTFtagExonParentTranscript} \
            --sjdbGTFtagExonParentGene ~{sjdbGTFtagExonParentGene} \
            --sjdbGTFtagExonParentGeneName ~{sjdbGTFtagExonParentGeneName} \
            --sjdbGTFtagExonParentGeneType ~{sjdbGTFtagExonParentGeneType} \
            --genomeChrBinNbits ~{genomeChrBinNbits} \
            --genomeSAindexNbases ~{genomeSAindexNbases} \
            --genomeSAsparseD ~{genomeSAsparseD} \
            --genomeSuffixLengthMax ~{genomeSuffixLengthMax} \
            --sjdbOverhang ~{sjdbOverhang} 

        rm "$gtf_name" "$ref_fasta"

        tar -czf ~{star_db_tar_gz} ~{db_name}
    >>>

    output {
        File star_db = star_db_tar_gz
    }

    runtime {
        cpu: ncpu
        memory: "~{memory_gb} GB"
        disk: "~{disk_size_gb} GB"
        docker: 'ghcr.io/stjudecloud/star:2.7.10a-0'
        maxRetries: max_retries
    }
}

task alignment {
    meta {
        description: "This WDL task runs the STAR aligner on a set of RNA-Seq FASTQ files."
    }

    parameter_meta {
        read_one_fastqs_gz: "An array of gzipped FASTQ files containing read one information"
        star_db_tar_gz: "A gzipped TAR file containing the STAR reference files. The name of the root directory which was archived must match the archive's filename without the `.tar.gz` extension."
        prefix: "Prefix for the BAM and other STAR files. The extensions `.Aligned.out.bam`, `.Log.final.out`, `.SJ.out.tab`, and `.Chimeric.out.junction` will be added."
        read_groups: "A string containing the read group information to output in the BAM file. If including multiple read group fields per-read group, they should be space delimited. Read groups should be comma separated, with a space on each side (i.e. ' , '). The ID field must come first for each read group and must match the basename of a fastq file (up to the first period). Example: `ID:rg1 PU:flowcell1.lane1 SM:sample1 PL:illumina LB:sample1_lib1 , ID:rg2 PU:flowcell1.lane2 SM:sample1 PL:illumina LB:sample1_lib1`"
        read_two_fastqs_gz: "An array of FASTQ files containing read two information"
        outSJfilterIntronMaxVsReadN: "maximum gap allowed for junctions supported by 1,2,3,,,N reads. i.e. by default junctions supported by 1 read can have gaps <=50000b, by 2 reads: <=100000b, by 3 reads: <=200000b. by >=4 reads any gap <=alignIntronMax. Does not apply to annotated junctions."
        outSJfilterOverhangMin: "minimum overhang length for splice junctions on both sides for: (1) non-canonical motifs, (2) GT/AG and CT/AC motif, (3) GC/AG and CT/GC motif, (4) AT/AC and GT/AT motif. -1 means no output for that motif. Does not apply to annotated junctions."
        outSJfilterCountUniqueMin: "minimum uniquely mapping read count per junction for: (1) non-canonical motifs, (2) GT/AG and CT/AC motif, (3) GC/AG and CT/GC motif, (4) AT/AC and GT/AT motif. -1 means no output for that motif. Junctions are output if one of outSJfilterCountUniqueMin *OR* outSJfilterCountTotalMin conditions are satisfied. Does not apply to annotated junctions."
        outSJfilterCountTotalMin: "minimum total (multi-mapping+unique) read count per junction for: (1) non-canonical motifs, (2) GT/AG and CT/AC motif, (3) GC/AG and CT/GC motif, (4) AT/AC and GT/AT motif. -1 means no output for that motif. Junctions are output if one of outSJfilterCountUniqueMin *OR* outSJfilterCountTotalMin conditions are satisfied. Does not apply to annotated junctions."
        outSJfilterDistToOtherSJmin: "minimum allowed distance to other junctions' donor/acceptor for: (1) non-canonical motifs, (2) GT/AG and CT/AC motif, (3) GC/AG and CT/GC motif, (4) AT/AC and GT/AT motif. Does not apply to annotated junctions."
        alignSJstitchMismatchNmax: "maximum number of mismatches for stitching of the splice junctions (-1: no limit) for: (1) non-canonical motifs, (2) GT/AG and CT/AC motif, (3) GC/AG and CT/GC motif, (4) AT/AC and GT/AT motif"
        clip3pAdapterSeq: {
            description: "adapter sequences to clip from 3p of each mate. `left` applies to read one and `right` applies to read two.",
            choices: {
                None: "No 3p adapter trimming will be performed",
                sequence: "A nucleotide sequence string of any length, matching the regex `/[ATCG]+/`",
                polyA: "polyA sequence with the length equal to read length"
            }
        }
        clip3pAdapterMMp: "max proportion of mismatches for 3p adapter clipping for each mate. `left` applies to read one and `right` applies to read two."
        alignEndsProtrude: {
            description: "allow protrusion of alignment ends, i.e. start (end) of the +strand mate downstream of the start (end) of the -strand mate. `left`: maximum number of protrusion bases allowed. `right`: see `choices` below.",
            choices: {
                ConcordantPair: "report alignments with non-zero protrusion as concordant pairs",
                DiscordantPair: "report alignments with non-zero protrusion as discordant pairs"
            }
        }
        clip3pNbases: "number of bases to clip from 3p of each mate. `left` applies to read one and `right` applies to read two."
        clip3pAfterAdapterNbases: "number of bases to clip from 3p of each mate after the adapter clipping. `left` applies to read one and `right` applies to read two."
        clip5pNbases: "number of bases to clip from 5p of each mate. `left` applies to read one and `right` applies to read two."
        readNameSeparator: "character(s) separating the part of the read names that will be trimmed in output (read name after space is always trimmed)"
        clipAdapterType: {
            description: "adapter clipping type",
            choices: {
                Hamming: "adapter clipping based on Hamming distance, with the number of mismatches controlled by --clip5pAdapterMMp",
                CellRanger4: "5p and 3p adapter clipping similar to CellRanger4. Utilizes Opal package by Martin Šošić: https://github.com/Martinsos/opal",
                None: "no adapter clipping, all other clip* parameters are disregarded"
            }
        }
        outSAMstrandField: {
            description: "Cufflinks-like strand field flag",
            choices: {
                None: "not used",
                intronMotif: "strand derived from the intron motif. This option changes the output alignments: reads with inconsistent and/or non-canonical introns are filtered out."
            }
        }
        outSAMattributes: {
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
                cN: "number of bases clipped from the read ends: 5' and 3'"
            }
        }
        outSAMunmapped: {
            description: "output of unmapped reads in the SAM format.",
            choices: {
                None: "no output **[STAR default]**",
                Within: "output unmapped reads within the main SAM file (i.e. Aligned.out.sam) **[WDL default]**"
            }
        }
        outSAMorder: {
            description: "type of sorting for the SAM output",
            choices: {
                Paired: "one mate after the other for all paired alignments",
                PairedKeepInputOrder: "one mate after the other for all paired alignments, the order is kept the same as in the input FASTQ files"
            }
        }
        outSAMreadID: {
            description: "read ID record type",
            choices: {
                Standard: "first word (until space) from the FASTx read ID line, removing /1,/2 from the end",
                Number: "read number (index) in the FASTx file"
            }
        }
        outSAMtlen: {
            description: "calculation method for the TLEN field in the SAM/BAM files",
            choices: {
                left_plus: "leftmost base of the (+)strand mate to rightmost base of the (-)mate. (+)sign for the (+)strand mate",
                left_any: "leftmost base of any mate to rightmost base of any mate. (+)sign for the mate with the leftmost base. This is different from `left_plus` for overlapping mates with protruding ends"
            }
        }
        outFilterType: {
            description: "type of filtering",
            choices: {
                Normal: "standard filtering using only current alignment",
                BySJout: "keep only those reads that contain junctions that passed filtering into SJ.out.tab"
            }
        }
        outFilterIntronMotifs: {
            description: "filter alignment using their motifs",
            choices: {
                None: "no filtering",
                RemoveNoncanonical: "filter out alignments that contain non-canonical junctions",
                RemoveNoncanonicalUnannotated: "filter out alignments that contain non-canonical unannotated junctions when using annotated splice junctions database. The annotated non-canonical junctions will be kept."
            }
        }
        outFilterIntronStrands: {
            description: "filter alignments",
            choices: {
                None: "no filtering",
                RemoveInconsistentStrands: "remove alignments that have junctions with inconsistent strands"
            }
        }
        outSJfilterReads: {
            description: "which reads to consider for collapsed splice junctions output",
            choices: {
                All: "all reads, unique- and multi-mappers",
                Unique: "uniquely mapping reads only"
            }
        }
        alignEndsType: {
            description: "type of read ends alignment",
            choices: {
                Local: "standard local alignment with soft-clipping allowed",
                EndToEnd: "force end-to-end read alignment, do not soft-clip",
                Extend5pOfRead1: "fully extend only the 5p of the read1, all other ends: local alignment",
                Extend5pOfReads12: "fully extend only the 5p of the both read1 and read2, all other ends: local alignment"
            }
        }
        alignSoftClipAtReferenceEnds: {
            description: "allow the soft-clipping of the alignments past the end of the chromosomes",
            choices: {
                Yes: "allow",
                No: "prohibit, useful for compatibility with Cufflinks"
            }
        }
        alignInsertionFlush: {
            description: "how to flush ambiguous insertion positions",
            choices: {
                None: "insertions are not flushed",
                Right: "insertions are flushed to the right"
            }
        }
        chimOutType: {
            description: "type of chimeric output",
            choices: {
                Junctions: "Chimeric.out.junction",
                WithinBAM_HardClip: "output into main aligned BAM files (Aligned.*.bam). Hard-clipping in the CIGAR for supplemental chimeric alignments."
                WithinBAM_SoftClip: "output into main aligned BAM files (Aligned.*.bam). Soft-clipping in the CIGAR for supplemental chimeric alignments."
            }
        }
        chimFilter: {
            description: "different filters for chimeric alignments",
            choices: {
                None: "no filtering",
                banGenomicN: "Ns are not allowed in the genome sequence around the chimeric junction"
            }
        }
        chimOutJunctionFormat: {
            description: "formatting type for the Chimeric.out.junction file",
            choices: {
                plain: "no comment lines/headers",
                comments: "comment lines at the end of the file: command line and Nreads: total, unique/multi-mapping"
            }
        }
        twopassMode: {
            description: "2-pass mapping mode",
            choices: {
                None: "1-pass mapping **[STAR default]**",
                Basic: "basic 2-pass mapping, with all 1st pass junctions inserted into the genome indices on the fly **[WDL default]**"
            }
        }
        use_all_cores: "Use all cores? Recommended for cloud environments. Not recommended for cluster environments."
        outFilterMismatchNoverLmax: "alignment will be output only if its ratio of mismatches to *mapped* length is less than or equal to this value"
        outFilterMismatchNoverReadLmax: "alignment will be output only if its ratio of mismatches to *read* length is less than or equal to this value"
        outFilterScoreMinOverLread: "same as outFilterScoreMin, but normalized to read length (sum of mates' lengths for paired-end reads)"
        outFilterMatchNminOverLread: "same as outFilterMatchNmin, but normalized to the read length (sum of mates' lengths for paired-end reads)"
        scoreGenomicLengthLog2scale: "extra score logarithmically scaled with genomic length of the alignment: scoreGenomicLengthLog2scale*log2(genomicLength)"
        seedSearchStartLmaxOverLread: "seedSearchStartLmax normalized to read length (sum of mates' lengths for paired-end reads)"
        alignSplicedMateMapLminOverLmate: "alignSplicedMateMapLmin normalized to mate length"
        peOverlapMMp: "maximum proportion of mismatched bases in the overlap area"
        runRNGseed: "random number generator seed"
        sjdbScore: "extra alignment score for alignments that cross database junctions"
        readMapNumber: "number of reads to map from the beginning of the file. -1 to map all reads"
        readQualityScoreBase: "number to be subtracted from the ASCII code to get Phred quality score"
        limitOutSJoneRead: "max number of junctions for one read (including all multi-mappers)"
        limitOutSJcollapsed: "max number of collapsed junctions"
        limitSjdbInsertNsj: "maximum number of junction to be inserted to the genome on the fly at the mapping stage, including those from annotations and those detected in the 1st step of the 2-pass run"
        outQSconversionAdd: "add this number to the quality score (e.g. to convert from Illumina to Sanger, use -31)"
        outSAMattrIHstart: "start value for the IH attribute. 0 may be required by some downstream software, such as Cufflinks or StringTie."
        outSAMmapqUnique: "`0-255`: the MAPQ value for unique mappers"
        outSAMflagOR: "`0-65535`: sam FLAG will be bitwise OR'd with this value, i.e. FLAG=FLAG | outSAMflagOR. This is applied after all flags have been set by STAR, and after outSAMflagAND. Can be used to set specific bits that are not set otherwise."
        outSAMflagAND: "`0-65535`: sam FLAG will be bitwise AND'd with this value, i.e. FLAG=FLAG & outSAMflagOR. This is applied after all flags have been set by STAR, but before outSAMflagOR. Can be used to unset specific bits that are not set otherwise."
        outFilterMultimapScoreRange: "the score range below the maximum score for multimapping alignments"
        outFilterMultimapNmax: "maximum number of loci the read is allowed to map to. Alignments (all of them) will be output only if the read maps to no more loci than this value. Otherwise no alignments will be output, and the read will be counted as 'mapped to too many loci' in the Log.final.out. **[STAR default]**: `10`. **[WDL default]**: `20`."
        outFilterMismatchNmax: "alignment will be output only if it has no more mismatches than this value"
        outFilterScoreMin: "alignment will be output only if its score is higher than or equal to this value"
        outFilterMatchNmin: "alignment will be output only if the number of matched bases is higher than or equal to this value"
        scoreGap: "splice junction penalty (independent on intron motif)"
        scoreGapNoncan: "non-canonical junction penalty (in addition to scoreGap)"
        scoreGapGCAG: "GC/AG and CT/GC junction penalty (in addition to scoreGap)"
        scoreGapATAC: "AT/AC and GT/AT junction penalty (in addition to scoreGap)"
        scoreDelOpen: "deletion open penalty"
        scoreDelBase: "deletion extension penalty per base (in addition to scoreDelOpen)"
        scoreInsOpen: "insertion open penalty"
        scoreInsBase: "insertion extension penalty per base (in addition to scoreInsOpen)"
        scoreStitchSJshift: "maximum score reduction while searching for SJ boundaries in the stitching step"
        seedSearchStartLmax: "defines the search start point through the read - the read is split into pieces no longer than this value"
        seedSearchLmax: "defines the maximum length of the seeds, if =0 seed length is not limited"
        seedMultimapNmax: "only pieces that map fewer than this value are utilized in the stitching procedure"
        seedPerReadNmax: "max number of seeds per read"
        seedPerWindowNmax: "max number of seeds per window"
        seedNoneLociPerWindow: "max number of one seed loci per window"
        seedSplitMin: "min length of the seed sequences split by Ns or mate gap"
        seedMapMin: "min length of seeds to be mapped"
        alignIntronMin: "minimum intron size: genomic gap is considered intron if its length>=alignIntronMin, otherwise it is considered Deletion"
        alignIntronMax: "maximum intron size, if 0, max intron size will be determined by (2^winBinNbits)*winAnchorDistNbins. **[STAR default]**: `0`. **[WDL default]**: `500000`."
        alignMatesGapMax: "maximum gap between two mates, if 0, max intron gap will be determined by (2^winBinNbits)*winAnchorDistNbins. **[STAR default]**: `0`. **[WDL default]**: `1000000`"
        alignSJoverhangMin: "minimum overhang (i.e. block size) for spliced alignments"
        alignSJDBoverhangMin: "minimum overhang (i.e. block size) for annotated (sjdb) spliced alignments. **[STAR default]**: `3`. **[WDL default]**: `1`."
        alignSplicedMateMapLmin: "minimum mapped length for a read mate that is spliced"
        alignWindowsPerReadNmax: "max number of windows per read"
        alignTranscriptsPerWindowNmax: "max number of transcripts per window"
        alignTranscriptsPerReadNmax: "max number of different alignments per read to consider"
        peOverlapNbasesMin: "minimum number of overlap bases to trigger mates merging and realignment. Specify >0 value to switch on the 'merging of overlapping mates' algorithm."
        winAnchorMultimapNmax: "max number of loci anchors are allowed to map to"
        winBinNbits: "=log2(winBin), where winBin is the size of the bin for the windows/clustering, each window will occupy an integer number of bins"
        winAnchorDistNbins: "max number of bins between two anchors that allows aggregation of anchors into one window"
        winFlankNbins: "=log2(winFlank), where winFlank is the size of the left and right flanking regions for each window"
        chimSegmentMin: "minimum length of chimeric segment length, if ==0, no chimeric output"
        chimScoreMin: "minimum total (summed) score of the chimeric segments"
        chimScoreDropMax: "max drop (difference) of chimeric score (the sum of scores of all chimeric segments) from the read length"
        chimScoreSeparation: "minimum difference (separation) between the best chimeric score and the next one"
        chimScoreJunctionNonGTAG: "penalty for a non-GT/AG chimeric junction"
        chimJunctionOverhangMin: "minimum overhang for a chimeric junction"
        chimSegmentReadGapMax: "maximum gap in the read sequence between chimeric segments"
        chimMainSegmentMultNmax: "maximum number of multi-alignments for the main chimeric segment. =1 will prohibit multimapping main segments."
        chimMultimapNmax: "maximum number of chimeric multi-alignments. `0`: use the old scheme for chimeric detection which only considered unique alignments"
        chimMultimapScoreRange: "the score range for multi-mapping chimeras below the best chimeric score. Only works with --chimMultimapNmax > 1."
        chimNonchimScoreDropMin: "to trigger chimeric detection, the drop in the best non-chimeric alignment score with respect to the read length has to be greater than this value"
        twopass1readsN: "number of reads to process for the 1st step. Use default (`-1`) to map all reads in the first step"
        ncpu: "Number of cores to allocate for task"
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
        max_retries: "Number of times to retry in case of failure"
    }

    input {
        File star_db_tar_gz
        Array[File] read_one_fastqs_gz
        String prefix
        String? read_groups
        Array[File] read_two_fastqs_gz = []
        Array[Int] outSJfilterIntronMaxVsReadN = [50000, 100000, 200000]
        SJ_motifs outSJfilterOverhangMin = SJ_motifs {
            noncanonical_motifs: 30,
            GT_AG_and_CT_AC_motif: 12,
            GC_AG_and_CT_GC_motif: 12,
            AT_AC_and_GT_AT_motif: 12
        }
        SJ_motifs outSJfilterCountUniqueMin = SJ_motifs {
            noncanonical_motifs: 3,
            GT_AG_and_CT_AC_motif: 1,
            GC_AG_and_CT_GC_motif: 1,
            AT_AC_and_GT_AT_motif: 1
        }
        SJ_motifs outSJfilterCountTotalMin = SJ_motifs {
            noncanonical_motifs: 3,
            GT_AG_and_CT_AC_motif: 1,
            GC_AG_and_CT_GC_motif: 1,
            AT_AC_and_GT_AT_motif: 1
        }
        SJ_motifs outSJfilterDistToOtherSJmin = SJ_motifs {
            noncanonical_motifs: 10,
            GT_AG_and_CT_AC_motif: 0,
            GC_AG_and_CT_GC_motif: 5,
            AT_AC_and_GT_AT_motif: 10
        }
        SJ_motifs alignSJstitchMismatchNmax = SJ_motifs {
            noncanonical_motifs: 0,
            GT_AG_and_CT_AC_motif: -1,
            GC_AG_and_CT_GC_motif: 0,
            AT_AC_and_GT_AT_motif: 0
        }
        Pair[String, String] clip3pAdapterSeq = ("None", "None")
        Pair[Float, Float] clip3pAdapterMMp = (0.1, 0.1)
        Pair[Int, String] alignEndsProtrude = (0, "ConcordantPair")
        Pair[Int, Int] clip3pNbases = (0, 0)
        Pair[Int, Int] clip3pAfterAdapterNbases = (0, 0)
        Pair[Int, Int] clip5pNbases = (0, 0)
        String readNameSeparator = "/"
        String clipAdapterType = "Hamming"
        String outSAMstrandField = "intronMotif"
        String outSAMattributes = "NH HI AS nM NM MD XS"
        String outSAMunmapped = "Within"
        String outSAMorder = "Paired"
        String outSAMreadID = "Standard"
        String outSAMtlen = "left_plus"
        String outFilterType = "Normal"
        String outFilterIntronMotifs = "None"
        String outFilterIntronStrands = "RemoveInconsistentStrands"
        String outSJfilterReads = "All"
        String alignEndsType = "Local"
        String alignSoftClipAtReferenceEnds = "Yes"
        String alignInsertionFlush = "None"
        String chimOutType = "Junctions"
        String chimFilter = "banGenomicN"
        String chimOutJunctionFormat = "plain"
        String twopassMode = "Basic"
        Boolean use_all_cores = false
        Float outFilterMismatchNoverLmax = 0.3
        Float outFilterMismatchNoverReadLmax = 1.0
        Float outFilterScoreMinOverLread = 0.66
        Float outFilterMatchNminOverLread = 0.66
        Float scoreGenomicLengthLog2scale = -0.25
        Float seedSearchStartLmaxOverLread = 1.0
        Float alignSplicedMateMapLminOverLmate = 0.66
        Float peOverlapMMp = 0.01
        Int runRNGseed = 777
        Int sjdbScore = 2
        Int readMapNumber = -1
        Int readQualityScoreBase = 33
        Int limitOutSJoneRead = 1000
        Int limitOutSJcollapsed = 1000000
        Int limitSjdbInsertNsj = 1000000
        Int outQSconversionAdd = 0
        Int outSAMattrIHstart = 1
        Int outSAMmapqUnique = 255
        Int outSAMflagOR = 0
        Int outSAMflagAND = 65535
        Int outFilterMultimapScoreRange = 1
        Int outFilterMultimapNmax = 20
        Int outFilterMismatchNmax = 10
        Int outFilterScoreMin = 0
        Int outFilterMatchNmin = 0
        Int scoreGap = 0
        Int scoreGapNoncan = -8
        Int scoreGapGCAG = -4
        Int scoreGapATAC = -8
        Int scoreDelOpen = -2
        Int scoreDelBase = -2
        Int scoreInsOpen = -2
        Int scoreInsBase = -2
        Int scoreStitchSJshift = 1
        Int seedSearchStartLmax = 50
        Int seedSearchLmax = 0
        Int seedMultimapNmax = 10000
        Int seedPerReadNmax = 1000
        Int seedPerWindowNmax = 50
        Int seedNoneLociPerWindow = 10
        Int seedSplitMin = 12
        Int seedMapMin = 5
        Int alignIntronMin = 21
        Int alignIntronMax = 500000
        Int alignMatesGapMax = 1000000
        Int alignSJoverhangMin = 5
        Int alignSJDBoverhangMin = 1
        Int alignSplicedMateMapLmin = 0
        Int alignWindowsPerReadNmax = 10000
        Int alignTranscriptsPerWindowNmax = 100
        Int alignTranscriptsPerReadNmax = 10000
        Int peOverlapNbasesMin = 0
        Int winAnchorMultimapNmax = 50
        Int winBinNbits = 16
        Int winAnchorDistNbins = 9
        Int winFlankNbins = 4
        Int chimSegmentMin = 0
        Int chimScoreMin = 0
        Int chimScoreDropMax = 20
        Int chimScoreSeparation = 10
        Int chimScoreJunctionNonGTAG = -1
        Int chimJunctionOverhangMin = 20
        Int chimSegmentReadGapMax = 0
        Int chimMainSegmentMultNmax = 10
        Int chimMultimapNmax = 0
        Int chimMultimapScoreRange = 1
        Int chimNonchimScoreDropMin = 20
        Int twopass1readsN = -1
        Int ncpu = 1
        Int memory_gb = 50
        Int modify_disk_size_gb = 0
        Int max_retries = 1
    }
    
    String star_db_dir = basename(star_db_tar_gz, ".tar.gz")

    Float read_one_fastqs_size = size(read_one_fastqs_gz, "GiB")
    Float read_two_fastqs_size = size(read_two_fastqs_gz, "GiB")
    Float star_db_tar_gz_size = size(star_db_tar_gz, "GiB")
    Int disk_size_gb = (
        (
            ceil(read_one_fastqs_size + read_two_fastqs_size + star_db_tar_gz_size) * 3
        ) + 10 + modify_disk_size_gb
    )

    Array[File] empty_array = []  # odd construction forced by WDL v1.1 spec

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi

        tar -xzf ~{star_db_tar_gz}

        # odd constructions a combination of needing white space properly parsed
        # and limitations of the WDL v1.1 spec
        python3 /home/sort_star_input.py \
            --read-one-fastqs "~{sep(',', read_one_fastqs_gz)}" \
            ~{if (read_two_fastqs_gz != empty_array) then "--read-two-fastqs" else ""} "~{
                sep(',', (
                    if (read_two_fastqs_gz != empty_array)
                    then read_two_fastqs_gz
                    else []
                ))
            }" \
            ~{if defined(read_groups) then "--read-groups" else ""} "~{
                if defined(read_groups)
                then read_groups
                else ''
            }"

        read -ra read_one_args < read_one_fastqs_sorted.txt
        read -ra read_two_args < read_two_fastqs_sorted.txt
        read -ra read_group_args < read_groups_sorted.txt
        STAR --readFilesIn "${read_one_args[@]}" "${read_two_args[@]}" \
            --readFilesCommand "gunzip -c" \
            --genomeDir ~{star_db_dir} \
            --runThreadN "$n_cores" \
            --outSAMtype BAM Unsorted \
            --outMultimapperOrder Random \
            --outFileNamePrefix ~{prefix + "."} \
            --twopassMode ~{twopassMode} \
            --outSAMattrRGline "${read_group_args[@]}" \
            --outSJfilterIntronMaxVsReadN ~{
                sep(' ', quote(outSJfilterIntronMaxVsReadN))
            } \
            --outSJfilterOverhangMin ~{sep(' ', quote([
                outSJfilterOverhangMin.noncanonical_motifs,
                outSJfilterOverhangMin.GT_AG_and_CT_AC_motif,
                outSJfilterOverhangMin.GC_AG_and_CT_GC_motif,
                outSJfilterOverhangMin.AT_AC_and_GT_AT_motif
            ]))} \
            --outSJfilterCountUniqueMin ~{sep(' ', quote([
                outSJfilterCountUniqueMin.noncanonical_motifs,
                outSJfilterCountUniqueMin.GT_AG_and_CT_AC_motif,
                outSJfilterCountUniqueMin.GC_AG_and_CT_GC_motif,
                outSJfilterCountUniqueMin.AT_AC_and_GT_AT_motif
            ]))} \
            --outSJfilterCountTotalMin ~{sep(' ', quote([
                outSJfilterCountTotalMin.noncanonical_motifs,
                outSJfilterCountTotalMin.GT_AG_and_CT_AC_motif,
                outSJfilterCountTotalMin.GC_AG_and_CT_GC_motif,
                outSJfilterCountTotalMin.AT_AC_and_GT_AT_motif
            ]))} \
            --outSJfilterDistToOtherSJmin ~{sep(' ', quote([
                outSJfilterDistToOtherSJmin.noncanonical_motifs,
                outSJfilterDistToOtherSJmin.GT_AG_and_CT_AC_motif,
                outSJfilterDistToOtherSJmin.GC_AG_and_CT_GC_motif,
                outSJfilterDistToOtherSJmin.AT_AC_and_GT_AT_motif
            ]))} \
            --alignSJstitchMismatchNmax ~{sep(' ', quote([
                alignSJstitchMismatchNmax.noncanonical_motifs,
                alignSJstitchMismatchNmax.GT_AG_and_CT_AC_motif,
                alignSJstitchMismatchNmax.GC_AG_and_CT_GC_motif,
                alignSJstitchMismatchNmax.AT_AC_and_GT_AT_motif
            ]))} \
            --clip3pAdapterSeq ~{clip3pAdapterSeq.left + ' ' + clip3pAdapterSeq.right} \
            --clip3pAdapterMMp ~{'~{clip3pAdapterMMp.left} ~{clip3pAdapterMMp.right}'} \
            --alignEndsProtrude ~{
                '~{alignEndsProtrude.left} ~{alignEndsProtrude.right}'
            } \
            --clip3pNbases ~{'~{clip3pNbases.left} ~{clip3pNbases.right}'} \
            --clip3pAfterAdapterNbases ~{
                '~{clip3pAfterAdapterNbases.left} ~{clip3pAfterAdapterNbases.right}'
            } \
            --clip5pNbases ~{'~{clip5pNbases.left} ~{clip5pNbases.right}'} \
            --readNameSeparator ~{readNameSeparator} \
            --clipAdapterType ~{clipAdapterType} \
            --outSAMstrandField ~{outSAMstrandField} \
            --outSAMattributes ~{outSAMattributes} \
            --outSAMunmapped ~{outSAMunmapped} \
            --outSAMorder ~{outSAMorder} \
            --outSAMreadID ~{outSAMreadID} \
            --outSAMtlen ~{
                if (outSAMtlen == "left_plus")
                then "1"
                else (
                    if (outSAMtlen == "left_any") then "2" else "error"
                )
            } \
            --outFilterType ~{outFilterType} \
            --outFilterIntronMotifs ~{outFilterIntronMotifs} \
            --outFilterIntronStrands ~{outFilterIntronStrands} \
            --outSJfilterReads ~{outSJfilterReads} \
            --alignEndsType ~{alignEndsType} \
            --alignSoftClipAtReferenceEnds ~{alignSoftClipAtReferenceEnds} \
            --alignInsertionFlush ~{alignInsertionFlush} \
            --chimOutType ~{chimOutType} \
            --chimFilter ~{chimFilter} \
            --chimOutJunctionFormat ~{chimOutJunctionFormat} \
            --outFilterMismatchNoverLmax ~{outFilterMismatchNoverLmax} \
            --outFilterMismatchNoverReadLmax ~{outFilterMismatchNoverReadLmax} \
            --outFilterScoreMinOverLread ~{outFilterScoreMinOverLread} \
            --outFilterMatchNminOverLread ~{outFilterMatchNminOverLread} \
            --scoreGenomicLengthLog2scale ~{scoreGenomicLengthLog2scale} \
            --seedSearchStartLmaxOverLread ~{seedSearchStartLmaxOverLread} \
            --alignSplicedMateMapLminOverLmate ~{alignSplicedMateMapLminOverLmate} \
            --peOverlapMMp ~{peOverlapMMp} \
            --runRNGseed ~{runRNGseed} \
            --sjdbScore ~{sjdbScore} \
            --readMapNumber ~{readMapNumber} \
            --readQualityScoreBase ~{readQualityScoreBase} \
            --limitOutSJoneRead ~{limitOutSJoneRead} \
            --limitOutSJcollapsed ~{limitOutSJcollapsed} \
            --limitSjdbInsertNsj ~{limitSjdbInsertNsj} \
            --outQSconversionAdd ~{outQSconversionAdd} \
            --outSAMattrIHstart ~{outSAMattrIHstart} \
            --outSAMmapqUnique ~{outSAMmapqUnique} \
            --outSAMflagOR ~{outSAMflagOR} \
            --outSAMflagAND ~{outSAMflagAND} \
            --outFilterMultimapScoreRange ~{outFilterMultimapScoreRange} \
            --outFilterMultimapNmax ~{outFilterMultimapNmax} \
            --outFilterMismatchNmax ~{outFilterMismatchNmax} \
            --outFilterScoreMin ~{outFilterScoreMin} \
            --outFilterMatchNmin ~{outFilterMatchNmin} \
            --scoreGap ~{scoreGap} \
            --scoreGapNoncan ~{scoreGapNoncan} \
            --scoreGapGCAG ~{scoreGapGCAG} \
            --scoreGapATAC ~{scoreGapATAC} \
            --scoreDelOpen ~{scoreDelOpen} \
            --scoreDelBase ~{scoreDelBase} \
            --scoreInsOpen ~{scoreInsOpen} \
            --scoreInsBase ~{scoreInsBase} \
            --scoreStitchSJshift ~{scoreStitchSJshift} \
            --seedSearchStartLmax ~{seedSearchStartLmax} \
            --seedSearchLmax ~{seedSearchLmax} \
            --seedMultimapNmax ~{seedMultimapNmax} \
            --seedPerReadNmax ~{seedPerReadNmax} \
            --seedPerWindowNmax ~{seedPerWindowNmax} \
            --seedNoneLociPerWindow ~{seedNoneLociPerWindow} \
            --seedSplitMin ~{seedSplitMin} \
            --seedMapMin ~{seedMapMin} \
            --alignIntronMin ~{alignIntronMin} \
            --alignIntronMax ~{alignIntronMax} \
            --alignMatesGapMax ~{alignMatesGapMax} \
            --alignSJoverhangMin ~{alignSJoverhangMin} \
            --alignSJDBoverhangMin ~{alignSJDBoverhangMin} \
            --alignSplicedMateMapLmin ~{alignSplicedMateMapLmin} \
            --alignWindowsPerReadNmax ~{alignWindowsPerReadNmax} \
            --alignTranscriptsPerWindowNmax ~{alignTranscriptsPerWindowNmax} \
            --alignTranscriptsPerReadNmax ~{alignTranscriptsPerReadNmax} \
            --peOverlapNbasesMin ~{peOverlapNbasesMin} \
            --winAnchorMultimapNmax ~{winAnchorMultimapNmax} \
            --winBinNbits ~{winBinNbits} \
            --winAnchorDistNbins ~{winAnchorDistNbins} \
            --winFlankNbins ~{winFlankNbins} \
            --chimSegmentMin ~{chimSegmentMin} \
            --chimScoreMin ~{chimScoreMin} \
            --chimScoreDropMax ~{chimScoreDropMax} \
            --chimScoreSeparation ~{chimScoreSeparation} \
            --chimScoreJunctionNonGTAG ~{chimScoreJunctionNonGTAG} \
            --chimJunctionOverhangMin ~{chimJunctionOverhangMin} \
            --chimSegmentReadGapMax ~{chimSegmentReadGapMax} \
            --chimMainSegmentMultNmax ~{chimMainSegmentMultNmax} \
            --chimMultimapNmax ~{chimMultimapNmax} \
            --chimMultimapScoreRange ~{chimMultimapScoreRange} \
            --chimNonchimScoreDropMin ~{chimNonchimScoreDropMin} \
            --twopass1readsN ~{twopass1readsN}
    >>>

    output {
        File star_log = prefix + ".Log.final.out"
        File star_bam = prefix + ".Aligned.out.bam"
        File star_junctions = prefix + ".SJ.out.tab"
        File? star_chimeric_junctions = prefix + ".Chimeric.out.junction"
    }

    runtime {
        cpu: ncpu
        memory: "~{memory_gb} GB"
        disk: "~{disk_size_gb} GB"
        docker: 'ghcr.io/stjudecloud/star:branch-star-2.7.10a-1'
        maxRetries: max_retries
    }
}

struct SJ_motifs {
    Int noncanonical_motifs
    Int GT_AG_and_CT_AC_motif
    Int GC_AG_and_CT_GC_motif
    Int AT_AC_and_GT_AT_motif
}
