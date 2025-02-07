#@ except: ContainerValue, CommentWhitespace

version 1.1

import "../data_structures/read_group.wdl"

task build {
    meta {
        description: "Builds a Bowtie2 index from a FASTA file"
        outputs: {
            index_files: "Array of Bowtie2 index files",
            bowtie_db_tar_gz: "Bowtie2 index files in tarred, gzipped format",
        }
    }

    parameter_meta {
        reference: "Reference genome in FASTA format (optionally gzipped)"
        prefix: "Prefix for the Bowtie2 index files"
        ncpu: "Number of threads to use"
    }

    input {
        File reference
        String prefix
        Int ncpu = 4
    }

    String base = basename(reference, ".gz")

    command <<<
        set -euo pipefail

        gunzip -c ~{reference} > ~{base} \
           || ln -sf ~{reference} ~{base}

        bowtie2-build --threads=~{ncpu} ~{base} ~{prefix}

        tar -czf ~{prefix}.tar.gz ~{prefix}.*
    >>>

    output {
        Array[File] index_files = glob("*.bt2")
        File bowtie_db_tar_gz = prefix + ".tar.gz"
    }

    runtime {
        cpu: ncpu
        memory: "20 GB"
        container: "ghcr.io/stjudecloud/bowtie2:branch-hic_workflow-2.5.4-0"
        maxRetries: 1
    }
}

# Several bowtie2 options are intentionally omitted from the task definition.
# For more information on these parameters see `meta.external_help`.
# These include:
# `--mm` as memory-mapped I/O does not work with WDL and containerization.
# `--no-hd` and `--no-sq` as they produce malformed BAM files.
# `--interleaved` as reads should be in separate files.
# `--sra-acc` as an internet connection cannot be guaranteed.
# `-b` because WDL does not handle exclusive options well.
# `-S` as output is passed through `samtools` and written as BAM.
# `--tab5` as it is an antiquated format.
# `--tab6` as it is an antiquated format.
# `--qseq` as it is an antiquated format.
# `-f` as reads cannot have quality scores in FASTA.
# `-f` as it is an antiquated format.
# `-c` as it is not recommended for use.
# `-F` as it is not conducive to quality BAM output.
# `--phred33` and `--phred64` as they are non-standard.
# `--solexa-quals` as it is non-standard.
# `--int-quals` as it is non-standard.
# `--very-fast`, `--fast`,
# `--sensitive`, and `--very-sensitive` as they are convenience options.
# `--very-fast-local`, `--fast-local`,
# `--sensitive-local`, and `--very-sensitive-local` as they are convenience options.
# `-n-ceil`
# `--align-paired-read` as BAM input is not supported.
# `--preserve-tags` as BAM input is not supported.
# `-o|--offrate` as this should be set when building the index.
# `--time` this is hardcoded to on so timing information is always available.
# `--quiet` as this is hardcoded to off so information is always captured in the logs.
# Other options only expose one option:
# `--un-gz` - omits `--un`, `--un-bz2`, and `--un-lz4`
# `--al-gz` - omits `--al`, `--al-bz2`, and `--al-lz4`
# `--un-conc-gz` - omits `--un-conc`, `--un-conc-bz2`, and `--un-conc-lz4`
# `--al-conc-gz` - omits `--al-conc`, `--al-conc-bz2`, and `--al-conc-lz4`
# `--trim-to` as `-5|--trim5` and `-3|--trim3` are supported instead.
# `--fr` - omits `--rf` and `--ff`
task align {
    meta {
        description: "Aligns reads to a reference genome using Bowtie2"
        external_help: "https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml"
        outputs: {
            aligned_bam: "Aligned reads in BAM format",
            unpaired_unaligned: "Unpaired reads that didn't align, in FASTQ format",
            unpaired_aligned: "Unpaired reads that aligned at least once, in FASTQ format",
            paired_discordant_read_one: "Read one from pairs that didn't align concordantly, in FASTQ format",
            paired_discordant_read_two: "Read two from pairs that didn't align concordantly, in FASTQ format",
            paired_concordant_read_one: "Read one from pairs that aligned concordantly at least once, in FASTQ format. Reads are also present in the BAM.",
            paired_concordant_read_two: "Read two from pairs that aligned concordantly at least once, in FASTQ format. Reads are also present in the BAM.",
            metrics: "Metrics file",
        }
    }

    parameter_meta {
        bowtie_db_tar_gz: "Bowtie2 index files in tarred, gzipped format"
        read_one_fastq_gz: "A gzipped FASTQ file containing read one information"
        read_two_fastq_gz: "A gzipped FASTQ file containing read two information"
        skip: "Skip the first N reads/read pairs"
        upto: "Only align the first N reads/read pairs"
        trim5: "Trim N bases from the 5' end of each read before alignment"
        trim3: "Trim N bases from the 3' end of each read before alignment"
        seed_mismatch: {
            description: "Maximum number of mismatches in the seed [0, 1] during multiseed alignment.",
            external_help: "https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#multiseed-heuristic",
            bowtie_option: "-N",
        }
        seed_substring: {
            description: "Length of seed substring. >3, <32",
            external_help: "https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#multiseed-heuristic",
            bowtie_option: "-L",
        }
        dpad: "Include N extra reference characters on sides of dynamic programming table. Set to 0 for gapless alignment."
        gbar: "Disallow gaps within N nucleotides of read ends"
        ignore_quals: "When calculating a mismatch penalty, always consider the quality value at the mismatched position to be the highest possible, regardless of the actual value. Treat all quality values as 30 on Phred scale."
        nofw: "Do not align to the forward (original) strand of the reference"
        norc: "Do not align to the reverse-complement strand of the reference"
        no_1mm_upfront: {
            description: "Do not allow 1 mismatch alignments before attempting to scan for the optimal seeded alignments",
            help: "By default, Bowtie 2 will attempt to find either an exact or a 1-mismatch end-to-end alignment for the read before trying the multiseed heuristic. Such alignments can be found very quickly, and many short read alignments have exact or near-exact end-to-end alignments. However, this can lead to unexpected alignments when the user also sets options governing the multiseed heuristic, like `seed_subtring` and `seed_mismatch`. For instance, if the user specifies `seed_mismatch` 0 and `seed_substring` equal to the length of the read, the user will be surprised to find 1-mismatch alignments reported. This option prevents Bowtie 2 from searching for 1-mismatch end-to-end alignments before using the multiseed heuristic, which leads to the expected behavior when combined with options such as `seed_substring` and `seed_mismatch`. This comes at the expense of speed.",
        }
        end_to_end: {
            description: "If true, entire read must align; no clipping. Else, local alignment; ends might be soft clipped",
            help: "If true, Bowtie 2 requires that the entire read align from one end to the other, without any trimming (or soft clipping) of characters from either end. The match bonus always equals 0 in this mode, so all alignment scores are less than or equal to 0, and the greatest possible alignment score is 0. If false, Bowtie 2 does not require that the entire read align from one end to the other. Rather, some characters may be omitted (soft clipped) from the ends in order to achieve the greatest possible alignment score. The match bonus is used in this mode, and the best possible alignment score is equal to the match bonus times the length of the read.",
        }
        match_bonus: {
            description: "bonus for match (0 for end-to-end or 2 for local)",
            bowtie2_option: "ma",
        }
        mismatch_penalty: {
            description: "max penalty for mismatch; lower quality = lower penalty",
            help: "Sets the maximum (MX) and minimum (MN) mismatch penalties, both integers. A number less than or equal to MX and greater than or equal to MN is subtracted from the alignment score for each position where a read character aligns to a reference character, the characters do not match, and neither is an N. If `ignore_quals` is specified, the number subtracted quals MX. Otherwise, the number subtracted is MN + floor( (MX-MN)(MIN(Q, 40.0)/40.0) ) where Q is the Phred quality value.",
            bowtie2_option: "mp",
        }
        non_actg_penalty: "penalty for non-A/C/G/Ts (e.g. ambiguous character `N`) in read/reference"
        read_gap_open_extend: "Sets the read gap open and extend penalties. A read gap of length N gets a penalty of <int1> + N * <int2>."
        ref_gap_open_extend: "Sets the reference gap open and extend penalties. A reference gap of length N gets a penalty of <int1> + N * <int2>."
        max_aln_report: {
            description: "Report up to N alignments per read; MAPQ not meaningful. Mutualyl exclusive with `report_all_alignments`.",
            help: "When specified, bowtie2 searches for at most N distinct, valid alignments for each read. The search terminates when it can't find more distinct valid alignments, or when it finds <int>, whichever happens first. All alignments found are reported in descending order by alignment score. The alignment score for a paired-end alignment equals the sum of the alignment scores of the individual mates. Each reported read or pair alignment beyond the first has the SAM 'secondary' bit (which equals 256) set in its FLAGS field. For reads that have more than N distinct, valid alignments, bowtie2 does not guarantee that the N alignments reported are the best possible in terms of alignment score.",
        }
        report_all_alignments: "Report all alignments; very slow, MAPQ not meaningful. Mutually exclusive with `max_aln_report`."
        min_fragment_len: {
            description: "The minimum fragment length for valid paired-end alignments.",
            bowtie2_option: "I|minins",
        }
        max_fragment_len: {
            description: "The maximum fragment length for valid paired-end alignments.",
            bowtie2_option: "X|maxins",
        }
        no_mixed: "By default, when bowtie2 cannot find a concordant or discordant alignment for a pair, it then tries to find alignments for the individual mates. This option disables that behavior."
        no_discordant: "By default, bowtie2 looks for discordant alignments if it cannot find any concordant alignments. A discordant alignment is an alignment where both mates align uniquely, but that does not satisfy the paired-end constraints. This option disables that behavior."
        dovetail: {
            description: "If the mates 'dovetail', that is if one mate alignment extends past the beginning of the other such that the wrong mate begins upstream, consider that to be concordant.",
            external_help: "https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#mates-can-overlap-contain-or-dovetail-each-other",
        }
        no_contain: {
            description: "If one mate alignment contains the other, consider that to be non-concordant.",
            external_help: "https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#mates-can-overlap-contain-or-dovetail-each-other",
        }
        no_overlap: {
            description: "If one mate alignment overlaps the other at all, consider that to be non-concordant.",
            external_help: "https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#mates-can-overlap-contain-or-dovetail-each-other",
        }
        write_unpaired_unaligned: "Write unpaired reads that didn't align to a file in gzipped FASTQ format. These reads correspond to the SAM records with the FLAGS 0x4 bit set and neither the 0x40 nor 0x80 bits set."
        write_unpaired_aligned: "Write unpaired reads that aligned at least once to a file in gzipped FASTQ format. These reads correspond to the SAM records with the FLAGS 0x4, 0x40, and 0x80 bits unset."
        write_paired_discordant: "Write pairs that didn't align concordantly to a file in gzipped FASTQ format. These reads correspond to the SAM records with the FLAGS 0x4 bit set and either the 0x40 or 0x80 bit set (depending on whether it's mate #1 or #2)."
        write_paired_concordant: "Write pairs that aligned concordantly at least once to a file in gzipped FASTQ format. These reads correspond to the SAM records with the FLAGS 0x4 bit unset and either the 0x40 or 0x80 bit set (depending on whether it's mate #1 or #2)."
        metrics_file: "Write bowtie2 alignment metrics to file. Useful for debugging."
        metrics_stderr: "Write bowtie2 alignment metrics to STDERR. Not mutually exclusive with `metrics_file`."
        metrics_interval: "Report internal counters & metrics every N seconds. Only used if `metrics_file` and/or `metrics_stderr` is specified."
        no_unal: "Suppress SAM records for unaligned reads."
        read_group: "Read group record to include in output SAM/BAM header"
        addl_rg_text: "add <text> ('label:value') to @RG line of SAM header"
        omit_sec_seq: "When printing secondary alignments, Bowtie 2 by default will write out the SEQ and QUAL strings. Specifying this option causes Bowtie 2 to print an asterisk in those fields instead."
        sam_no_qname_trunc: "Suppress standard behavior of truncating readname at first whitespace at the expense of generating non-standard SAM."
        xeq: "Use '='/'X', instead of 'M,' to specify matches/mismatches in SAM record"
        soft_clipped_unmapped_tlen: "Exclude soft-clipped bases when reporting TLEN (template length). Only used if `end_to_end` is false."
        sam_append_comment: "Append FASTA/FASTQ comment to SAM record, where a comment is everything after the first space in the read name."
        sam_opt_config: "Toggle optional SAM fields. Prefix fields with `-` to turn off. Example: '-MD,YP,-AS' will disable the `MD` and `AS` fields and enable the `YP` field."
        ncpu: "Number of threads to use for alignment. Threads will synchronize when parsing reads and outputting alignments. Searching for alignments is highly parallel, and speedup is close to linear. Increasing threads increases Bowtie 2's memory footprint."
        reorder: {
            description: "Force SAM output order to match order of input reads",
            help: "Guarantees that output SAM records are printed in an order corresponding to the order of the reads in the original input file, even when `ncpu` is set greater than 1. Specifying `reorder` and setting `ncpu` greater than 1 causes Bowtie 2 to run somewhat slower and use somewhat more memory than if `reorder` were not specified. Has no effect if `ncpu` is set to 1, since output order will naturally correspond to input order in that case.",
        }
        qc_filter: {
            description: "Filter out reads that are bad according to QSEQ filter. Only has an effect when read format is --qseq.",
            common: false,
        }
        seed: "Seed for random number generator"
        non_deterministic: {
            description: "Seed random number generator arbitrarily instead of using read attributes",
            help: "Normally, Bowtie 2 re-initializes its pseudo-random generator for each read. It seeds the generator with a number derived from (a) the read name, (b) the nucleotide sequence, (c) the quality sequence, (d) the value of the `seed` option. This means that if two reads are identical (same name, same nucleotides, same qualities) Bowtie 2 will find and report the same alignment(s) for both, even if there was ambiguity. When `non_deterministic` is specified, Bowtie 2 re-initializes its pseudo-random generator for each read using the current time. This means that Bowtie 2 will not necessarily report the same alignment for two identical reads. This is counter-intuitive for some users, but might be more appropriate in situations where the input consists of many identical reads.",
        }
        prefix: "Prefix to use for output files"
        score_min: {
            description: "Sets a function governing the minimum alignment score needed for an alignment to be considered 'valid' (i.e. good enough to report).",
            external_help: "https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#setting-function-options",
        }
        interval_seed_substrings: {
            description: "Interval between seed substrings with respect to read len (S,1,1.15)",
            help: "Sets a function governing the interval between seed substrings to use during multiseed alignment. Since it's best to use longer intervals for longer reads, this parameter sets the interval as a function of the read length, rather than a single one-size-fits-all number.",
        }
        max_failed_extends: {
            description: "Give up extending after N failed extends in a row",
            help: "Up to N consecutive seed extension attempts can 'fail' before Bowtie 2 moves on, using the alignments found so far. A seed extension 'fails' if it does not yield a new best or a new second-best alignment. This limit is automatically adjusted up when `max_aln_report` or `report_all_alignments` are specified.",
        }
        repetitive_seeds: {
            description: "For reads with repetitive seeds, try N sets of seeds",
            help: "N is the maximum number of times Bowtie 2 will re-seed reads with repetitive seeds. When re-seeding, Bowtie 2 simply chooses a new set of reads (same length, same number of mismatches allowed) at different offsets and searches for more alignments. A read is considered to have repetitive seeds if the total number of seed hits divided by the number of seeds that aligned at least once is greater than 300.",
        }
        n_ceil: {
            description: "func for maximum number of non-A/C/G/Ts permitted in alignment",
            help: "Sets a function governing the maximum number of ambiguous characters (usually Ns and/or .s) allowed in a read as a function of read length. For instance, specifying -L,0,0.15 sets the N-ceiling function f to f(x) = 0 + 0.15 * x, where x is the read length. Reads exceeding this ceiling are filtered out."
        }
    }

    input {
        File bowtie_db_tar_gz
        File read_one_fastq_gz
        ReadGroup read_group
        String prefix
        File? read_two_fastq_gz
        String? addl_rg_text
        String? sam_opt_config
        Int? skip
        Int? upto
        Int? match_bonus
        Int? max_aln_report
        Bowtie2Function score_min = Bowtie2Function {
            function_type: "L",
            constant: -0.6,
            coefficient: -0.6,
        }
        Bowtie2Function interval_seed_substrings = Bowtie2Function {
            function_type: "S",
            constant: 1,
            coefficient: 1.15,
        }
        Bowtie2Function n_ceil = Bowtie2Function {
            function_type: "L",
            constant: 0,
            coefficient: 0.15,
        }
        Pair[Int, Int] read_gap_open_extend = (5, 3)
        Pair[Int, Int] ref_gap_open_extend = (5, 3)
        Pair[Int, Int] mismatch_penalty = (6, 2)
        Boolean end_to_end = true  # false: --local
        Boolean non_deterministic = false
        Boolean ignore_quals = false
        Boolean nofw = false
        Boolean norc = false
        Boolean no_1mm_upfront = false
        Boolean report_all_alignments = false
        Boolean no_mixed = false
        Boolean no_discordant = false
        Boolean dovetail = false
        Boolean no_contain = false
        Boolean no_overlap = false
        Boolean write_unpaired_unaligned = false
        Boolean write_unpaired_aligned = false
        Boolean write_paired_discordant = false
        Boolean write_paired_concordant = false
        Boolean metrics_file = false
        Boolean metrics_stderr = false
        Boolean no_unal = false
        Boolean omit_sec_seq = false
        Boolean sam_no_qname_trunc = false
        Boolean xeq = false
        Boolean soft_clipped_unmapped_tlen = false
        Boolean sam_append_comment = false
        Boolean reorder = false
        Boolean qc_filter = false
        Int trim5 = 0
        Int trim3 = 0
        Int seed_mismatch = 0
        Int seed_substring = 22
        Int dpad = 15
        Int gbar = 4
        Int non_actg_penalty = 1
        Int min_fragment_len = 0
        Int max_fragment_len = 500
        Int metrics_interval = 1
        Int ncpu = 8
        Int seed = 0
        Int max_failed_extends = 15
        Int repetitive_seeds = 2
    }

    #@ except: LineWidth
    command <<<
        set -euo pipefail

        # Setup reference files
        mkdir bowtie_db
        tar -C bowtie_db -xzf ~{bowtie_db_tar_gz} --no-same-owner
        PREFIX=$(basename bowtie_db/*.rev.1.bt2 ".rev.1.bt2")

        # Run bowtie2 alignment
        bowtie2 \
            ~{if defined(skip) then "--skip ~{skip}" else ""} \
            ~{if defined(upto) then "--upto ~{upto}" else ""} \
            --trim5 ~{trim5} \
            --trim3 ~{trim3} \
            -N ~{seed_mismatch} \
            -L ~{seed_substring} \
            --dpad ~{dpad} \
            --gbar ~{gbar} \
            ~{if ignore_quals then "--ignore-quals" else ""} \
            ~{if nofw then "--nofw" else ""} \
            ~{if norc then "--norc" else ""} \
            ~{if no_1mm_upfront then "--no-1mm-upfront" else ""} \
            ~{if end_to_end then "--end-to-end" else "--local"} \
            ~{if defined(match_bonus) then "--ma ~{match_bonus}" else ""} \
            --mp ~{mismatch_penalty.left},~{mismatch_penalty.right} \
            --np ~{non_actg_penalty} \
            --rdg ~{read_gap_open_extend.left},~{read_gap_open_extend.right} \
            --rfg ~{ref_gap_open_extend.left},~{ref_gap_open_extend.right} \
            ~{if defined(max_aln_report) then "-k ~{max_aln_report}" else ""} \
            ~{if report_all_alignments then "-a" else ""} \
            --minins ~{min_fragment_len} \
            --maxins ~{max_fragment_len} \
            ~{if no_mixed then "--no-mixed" else ""} \
            ~{if no_discordant then "--no-discordant" else ""} \
            ~{if dovetail then "--dovetail" else ""} \
            ~{if no_contain then "--no-contain" else ""} \
            ~{if no_overlap then "--no-overlap" else ""} \
            ~{(
                if write_unpaired_unaligned
                then "--un-gz ~{prefix}.unpaired_unaligned.fastq.gz"
                else ""
            )} \
            ~{(
                if write_unpaired_aligned
                then "--al-gz ~{prefix}.unpaired_aligned.fastq.gz"
                else ""
            )} \
            ~{(
                if write_paired_discordant
                then "--un-conc-gz ~{prefix}.paired_discordant.fastq.gz"
                else ""
            )} \
            ~{(
                if write_paired_concordant
                then "--al-conc-gz ~{prefix}.paired_concordant.fastq.gz"
                else ""
            )} \
            ~{if metrics_file then "--met-file ~{prefix}.metrics.txt" else ""} \
            ~{if metrics_stderr then "--met-stderr" else ""} \
            --met ~{metrics_interval} \
            ~{if no_unal then "--no-unal" else ""} \
            --rg-id ~{read_group.ID} \
            ~{if defined(read_group.BC) then "--rg BC:~{read_group.BC}" else ""} \
            ~{if defined(read_group.CN) then "--rg CN:~{read_group.CN}" else ""} \
            ~{if defined(read_group.DS) then "--rg DS:~{read_group.DS}" else ""} \
            ~{if defined(read_group.DT) then "--rg DT:~{read_group.DT}" else ""} \
            ~{if defined(read_group.FO) then "--rg FO:~{read_group.FO}" else ""} \
            ~{if defined(read_group.KS) then "--rg KS:~{read_group.KS}" else ""} \
            ~{if defined(read_group.LB) then "--rg LB:~{read_group.LB}" else ""} \
            ~{if defined(read_group.PG) then "--rg PG:~{read_group.PG}" else ""} \
            ~{if defined(read_group.PI) then "--rg PI:~{read_group.PI}" else ""} \
            ~{if defined(read_group.PL) then "--rg PL:~{read_group.PL}" else ""} \
            ~{if defined(read_group.PM) then "--rg PM:~{read_group.PM}" else ""} \
            ~{if defined(read_group.PU) then "--rg PU:~{read_group.PU}" else ""} \
            ~{if defined(read_group.SM) then "--rg SM:~{read_group.SM}" else ""} \
            ~{if defined(addl_rg_text) then "--rg ~{addl_rg_text}" else ""} \
            ~{if omit_sec_seq then "--omit-sec-seq" else ""} \
            ~{if sam_no_qname_trunc then "--sam-no-qname-trunc" else ""} \
            ~{if xeq then "--xeq" else ""} \
            ~{if soft_clipped_unmapped_tlen then "--soft-clipped-unmapped-tlen" else ""} \
            ~{if sam_append_comment then "--sam-append-comment" else ""} \
            ~{if defined(sam_opt_config) then "--sam-opt ~{sam_opt_config}" else ""} \
            -p ~{ncpu} \
            ~{if reorder then "--reorder" else ""} \
            ~{if qc_filter then "--qc-filter" else ""} \
            --seed ~{seed} \
            --time \
            ~{if non_deterministic then "--non-deterministic" else ""} \
            ~{"--score-min ~{score_min.function_type},~{score_min.constant},~{score_min.coefficient}"} \
            ~{"-i ~{interval_seed_substrings.function_type},~{interval_seed_substrings.constant},~{interval_seed_substrings.coefficient}"} \
            ~{"--n-ceil ~{n_ceil.function_type},~{n_ceil.constant},~{n_ceil.coefficient}"} \
            -D ~{max_failed_extends} \
            -R ~{repetitive_seeds} \
            -x bowtie_db/"$PREFIX" \
            ~{(
                if defined(read_two_fastq_gz)
                then "-1 ~{read_one_fastq_gz}"
                else "-U ~{read_one_fastq_gz}"
            )} \
            ~{if defined(read_two_fastq_gz) then "-2 ~{read_two_fastq_gz}" else ""} \
            | samtools view -bS - > ~{prefix}.bam    
    >>>

    output {
        File aligned_bam = prefix + ".bam"
        File? unpaired_unaligned = prefix + ".unpaired_unaligned.fastq.gz"
        File? unpaired_aligned = prefix + ".unpaired_aligned.fastq.gz"
        File? paired_discordant_read_one = prefix + ".paired_discordant.fastq.1.gz"
        File? paired_discordant_read_two = prefix + ".paired_discordant.fastq.2.gz"
        File? paired_concordant_read_one = prefix + ".paired_concordant.fastq.1.gz"
        File? paired_concordant_read_two = prefix + ".paired_concordant.fastq.2.gz"
        File? metrics = prefix + ".metrics.txt"
    }

    runtime {
        container: "ghcr.io/stjudecloud/bowtie2:branch-hic_workflow-2.5.4-0"
        cpu: ncpu
        memory: "20 GB"
    }
}

# Bowtie2 accepts several function parameters.
# The first term is a function type. Available function types are:
# - C - constant
# - L - linear
# - S - square-root
# - G - natural logarithm
# The constant and coefficient types may be negative and/or floating point numbers.
#
# An example input JSON entry might look like:
# ```
# {
#     "interval_seed_substrings": {
#         "function_type": "S",
#         "constant": 1,
#         "coefficient": 0.50,
#     },
# }
# ```
# See the function documentation in bowtie2
# (https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#setting-function-options)
struct Bowtie2Function {
    String function_type
    Float constant
    Float coefficient
}
