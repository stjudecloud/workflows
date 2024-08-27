version 1.1

import "../data_structures/read_group.wdl"

task build {
    meta {
        description: "Builds a Bowtie2 index from a FASTA file"
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
        #container: "ghcr.io/stjudecloud/bowtie2:2.5.4"
        container: "adthrasher/bowtie2:2.5.4"
        maxRetries: 1
    }
}

task align {
    meta {

    }

    parameter_meta {
        bowtie_db_tar_gz: "Bowtie2 index files in tarred, gzipped format"
        read_one_fastq_gz: "An array of gzipped FASTQ files containing read one information"
        read_two_fastq_gz: "An array of gzipped FASTQ files containing read two information"
        skip: "Skip the first N reads/read pairs"
        upto: "Only align the first N reads/read pairs"
        trim5: "Trim N bases from the 5' end of each read"
        trim3: "Trim N bases from the 3' end of each read"
        seed_mismatch: "Maximum number of mismatches in the seed [0, 1]"
        seed_substring: "Length of seed substring. >3, <32"
        dpad: "include N extra ref chars on sides of DP table"
        gbar: "disallow gaps within N nucs of read extremes"
        ignore_quals: "treat all quality values as 30 on Phred scale"
        nofw: "do not align forward (original) version of read"
        norc: "do not align reverse-complement version of read"
        no_1mm_upfront: "do not allow 1 mismatch alignments before attempting to scan for the optimal seeded alignments"
        end_to_end: "If true, entire read must align; no clipping. Else, local alignment; ends might be soft clipped"
        match_bonus: "bonus for match (0 for end-to-end or 2 for local)"
        mismatch_penalty: "max penalty for mismatch; lower qual = lower penalty"
        non_actg_penalty: "penalty for non-A/C/G/Ts in read/ref"
        read_gap_open_extend: "read gap open, extend penalties"
        ref_gap_open_extend: "ref gap open, extend penalties"
        max_aln_report: "report up to N alns per read; MAPQ not meaningful"
        report_all_alignments: "report all alignments; very slow, MAPQ not meaningful"
        min_fragment_len: "minimum fragment length"
        max_fragment_len: "maximum fragment length"
        no_mixed: "suppress unpaired alignments for paired reads"
        no_discordant: "suppress discordant alignments for paired reads"
        dovetail: "concordant when mates extend past each other"
        no_contain: "not concordant when one mate alignment contains other"
        no_overlap: "not concordant when mates overlap at all"
        time: "print wall-clock time taken by search phases"
        write_unpaired_unaligned: "write unpaired reads that didn't align"
        write_unpaired_aligned: "write unpaired reads that aligned at least once"
        write_paired_discordant: "write pairs that didn't align concordantly"
        write_paired_concordant: "write pairs that aligned concordantly at least once"
        quiet: "print nothing to stderr except serious errors"
        metrics_file: "write metrics to file"
        metrics_stderr: "write metrics to stderr"
        metrics_interval: "report internal counters & metrics every N seconds"
        no_unal: "suppress SAM records for unaligned reads"
        no_head: "suppress header lines, i.e. lines starting with @"
        no_sq: "suppress @SQ header lines"
        rg: "Read group record for SAM/BAM header"
        addl_rg_text: "add <text> ('label:value') to @RG line of SAM header"
        omit_sec_seq: "put '*' in SEQ and QUAL fields for secondary alignments."
        sam_no_quane_trunc: "Suppress standard behavior of truncating readname at first whitespace at the expense of generating non-standard SAM."
        xeq: "Use '='/'X', instead of 'M,' to specify matches/mismatches in SAM record"
        soft_clipped_unmapped_tlen: "Exclude soft-clipped bases when reporting TLEN."
        sam_append_comment: "Append FASTA/FASTQ comment to SAM record."
        sam_opt_config: "Use <config>, example '-MD,YP,-AS', to toggle SAM optional fields."
        threads: "Number of alignment threads to use"
        reorder: "force SAM output order to match order of input reads"
        memory_map: "use memory-mapped I/O for index; many 'bowtie's can share"
        qc_filter: "filter out reads that are bad according to QSEQ filter"
        seed: "seed for random number generator"
        non_deterministic: "seed random number generator arbitrarily instead of using read attributes"
        prefix: "Prefix to use for output files"
        score_min: "min acceptable alignment score w/r/t read length (G,20,8 for local, L,-0.6,-0.6 for end-to-end)"
        interval_seed_substrings: "interval between seed substrings w/r/t read len (S,1,1.15)"
        max_failed_extends: "give up extending after N failed extends in a row"
        repetitive_seeds: "for reads w/ repetitive seeds, try N sets of seeds"
    }

    input {
        File bowtie_db_tar_gz
        File read_one_fastq_gz
        File? read_two_fastq_gz
        Int? skip
        Int? upto
        Int trim5 = 0
        Int trim3 = 0
        Int seed_mismatch = 0
        Int seed_substring = 22
        Int dpad = 15
        Int gbar = 4
        Boolean ignore_quals = false
        Boolean nofw = false
        Boolean norc = false
        Boolean no_1mm_upfront = false
        Boolean end_to_end = true  # false: --local
        Int? match_bonus
        Int mismatch_penalty = 6
        Int non_actg_penalty = 1
        Pair[Int, Int] read_gap_open_extend = (5, 3)
        Pair[Int, Int] ref_gap_open_extend = (5, 3)
        Int? max_aln_report
        Boolean report_all_alignments = false
        Int min_fragment_len = 0
        Int max_fragment_len = 500
        # TODO: --fr, --rf, --ff args
        Boolean no_mixed = false
        Boolean no_discordant = false
        Boolean dovetail = false
        Boolean no_contain = false
        Boolean no_overlap = false
        Boolean time = false
        Boolean write_unpaired_unaligned = false
        Boolean write_unpaired_aligned = false
        Boolean write_paired_discordant = false
        Boolean write_paired_concordant = false
        Boolean quiet = false
        Boolean metrics_file = false
        Boolean metrics_stderr = false
        Int metrics_interval = 1
        Boolean no_unal = false
        Boolean no_head = false
        Boolean no_sq = false
        ReadGroup rg
        String? addl_rg_text
        Boolean omit_sec_seq = false
        Boolean sam_no_quane_trunc = false
        Boolean xeq = false
        Boolean soft_clipped_unmapped_tlen = false
        Boolean sam_append_comment = false
        String? sam_opt_config
        Int threads = 1
        Boolean reorder = false
        Boolean memory_map = false
        Boolean qc_filter = false
        Int seed = 0
        Boolean non_deterministic = false
        String prefix
        String? score_min
        String? interval_seed_substrings
        Int max_failed_extends = 15
        Int repetitive_seeds = 2
    }

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
            --mp ~{mismatch_penalty} \
            --np ~{non_actg_penalty} \
            --rdg ~{read_gap_open_extend.left},~{read_gap_open_extend.right} \
            --rfg ~{ref_gap_open_extend.left},~{ref_gap_open_extend.right} \
            ~{if defined(max_aln_report) then "-k ~{max_aln_report}" else ""} \
            ~{if report_all_alignments then "--all" else ""} \
            --minins ~{min_fragment_len} \
            --maxins ~{max_fragment_len} \
            ~{if no_mixed then "--no-mixed" else ""} \
            ~{if no_discordant then "--no-discordant" else ""} \
            ~{if dovetail then "--dovetail" else ""} \
            ~{if no_contain then "--no-contain" else ""} \
            ~{if no_overlap then "--no-overlap" else ""} \
            ~{if time then "--time" else ""} \
            ~{if write_unpaired_unaligned then "--un-gz ~{prefix}.unpaired_unaligned.fastq.gz" else ""} \
            ~{if write_unpaired_aligned then "--al-gz ~{prefix}.unpaired_aligned.gz" else ""} \
            ~{if write_paired_discordant then "--un-conc-gz ~{prefix}.paired_discordant.gz" else ""} \
            ~{if write_paired_concordant then "--al-conc-gz ~{prefix}.paired_concordant.gz" else ""} \
            ~{if quiet then "--quiet" else ""} \
            ~{if metrics_file then "--met-file ~{prefix}.metrics.txt" else ""} \
            ~{if metrics_stderr then "--met-stderr" else ""} \
            --met ~{metrics_interval} \
            ~{if no_unal then "--no-unal" else ""} \
            ~{if no_head then "--no-head" else ""} \
            ~{if no_sq then "--no-sq" else ""} \
            --rg-id ~{rg.ID} \
            ~{if defined(rg.BC) then "--rg BC:~{rg.BC}" else ""} \
            ~{if defined(rg.CN) then "--rg CN:~{rg.CN}" else ""} \
            ~{if defined(rg.DS) then "--rg DS:~{rg.DS}" else ""} \
            ~{if defined(rg.DT) then "--rg DT:~{rg.DT}" else ""} \
            ~{if defined(rg.FO) then "--rg FO:~{rg.FO}" else ""} \
            ~{if defined(rg.KS) then "--rg KS:~{rg.KS}" else ""} \
            ~{if defined(rg.LB) then "--rg LB:~{rg.LB}" else ""} \
            ~{if defined(rg.PG) then "--rg PG:~{rg.PG}" else ""} \
            ~{if defined(rg.PI) then "--rg PI:~{rg.PI}" else ""} \
            ~{if defined(rg.PL) then "--rg PL:~{rg.PL}" else ""} \
            ~{if defined(rg.PM) then "--rg PM:~{rg.PM}" else ""} \
            ~{if defined(rg.PU) then "--rg PU:~{rg.PU}" else ""} \
            ~{if defined(rg.SM) then "--rg SM:~{rg.SM}" else ""} \
            ~{if defined(addl_rg_text) then "--rg ~{addl_rg_text}" else ""} \
            ~{if omit_sec_seq then "--omit-sec-seq" else ""} \
            ~{if sam_no_quane_trunc then "--sam-no-qname-trunc" else ""} \
            ~{if xeq then "--xeq" else ""} \
            ~{if soft_clipped_unmapped_tlen then "--soft-clipped-unmapped-tlen" else ""} \
            ~{if sam_append_comment then "--sam-append-comment" else ""} \
            ~{if defined(sam_opt_config) then "--sam-opt ~{sam_opt_config}" else ""} \
            -p ~{threads} \
            ~{if reorder then "--reorder" else ""} \
            ~{if memory_map then "--mm" else ""} \
            ~{if qc_filter then "--qc-filter" else ""} \
            --seed ~{seed} \
            ~{if non_deterministic then "--non-deterministic" else ""} \
            ~{if defined(score_min) then "--score-min ~{score_min}" else ""} \
            ~{if defined(interval_seed_substrings) then "-i ~{interval_seed_substrings}" else ""} \
            -D ~{max_failed_extends} \
            -R ~{repetitive_seeds} \
            -x bowtie_db/"$PREFIX" \
            ~{if defined(read_two_fastq_gz) then "-1 ~{read_one_fastq_gz}" else "-U ~{read_one_fastq_gz}"} \
            ~{if defined(read_two_fastq_gz) then "-2 ~{read_two_fastq_gz}" else ""} \
            | samtools view -bS - > ~{prefix}.bam    
    >>>

    output {
        File aligned_bam = prefix + ".bam"
        File? unpaired_unaligned = prefix + ".unpaired_unaligned.fastq.gz"
        File? unpaired_aligned = prefix + ".unpaired_aligned.gz"
        File? paired_discordant = prefix + ".paired_discordant.gz"
        File? paired_concordant = prefix + ".paired_concordant.gz"
        File? metrics = prefix + ".metrics.txt"
    }

    runtime {
        #container: "ghcr.io/stjudecloud/bowtie2:2.5.4"
        container: "adthrasher/bowtie2:2.5.4"
        cpu: threads
        memory: "20 GB"
    }
}
