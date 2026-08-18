version 1.1

task build_salmon_index {
    meta {
        description: "Builds a Salmon index from a transcriptome FASTA file, for use in quantification"
        outputs: {
            salmon_index_tar_gz: "A gzipped TAR file containing the Salmon index files."
        }
    }

    parameter_meta {
        transcripts_fasta: "FASTA format file containing the reference transcriptome to index"
        index_name: {
            description: "Name for the output index, in compressed archive format. The suffix `.tar.gz` will be added.",
            group: "Common",
        }
        use_all_cores: {
            description: "Use all cores? Recommended for cloud environments.",
            group: "Resources",
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
        File transcripts_fasta
        String index_name = "salmon_index"
        Boolean use_all_cores = false
        Int ncpu = 4
        Int modify_disk_size_gb = 0
    }

    String salmon_index_filename = index_name + ".tar.gz"

    Float transcripts_fasta_size = size(transcripts_fasta, "GB")
    Int disk_size_gb = ceil(transcripts_fasta_size * 4) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi

        salmon index \
            -t "~{transcripts_fasta}" \
            -i "~{index_name}" \
            -p "$n_cores"

        tar -czf "~{salmon_index_filename}" "~{index_name}"
    >>>

    output {
        File salmon_index_tar_gz = salmon_index_filename
    }

    runtime {
        cpu: ncpu
        memory: "8 GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/salmon:1.9.0--h7e5ed60_0"
        maxRetries: 1
    }
}

task quant {
    meta {
        description: "Runs Salmon quant in mapping-based mode to quantify transcript-level expression from RNA-Seq reads, using a pre-built Salmon index"
        outputs: {
            quant_results_tar_gz: "A gzipped TAR file containing the Salmon quantification output directory, including `quant.sf`."
        }
    }

    parameter_meta {
        salmon_index_tar_gz: "A gzipped TAR file containing the Salmon index files. Suitable as the output of the `build_salmon_index` task."
        read_one_fastqs_gz: "An array of gzipped FASTQ files containing read one information"
        read_two_fastqs_gz: {
            description: "An array of gzipped FASTQ files containing read two information. Omit for single-end reads.",
            group: "Common",
        }
        lib_type: {
            description: "Salmon library type describing the relative orientation and strandedness of paired reads.",
            help: "Use `A` to let Salmon auto-detect the library type â€” recommended for most users.",
            group: "Common",
        }
        prefix: {
            description: "Prefix for the Salmon quantification output. The extension `.tar.gz` will be added.",
            group: "Common",
        }
        validate_mappings: {
            description: "Validate mappings using an alignment-based verification step.",
            group: "Salmon Options",
        }
        num_bootstraps: {
            description: "Salmon has the ability to optionally compute bootstrapped abundance estimates.",
            help: "This is done by resampling (with replacement) from the counts assigned to the fragment equivalence classes, and then re-running the optimization procedure for each such sample.",
            group: "Salmon Options",
        }
        incompat_prior: {
            description: "This parameter governs the a priori probability that a fragment mapping is nonetheless the correct mapping.",
            help: "Specifically, this is for a fragment mapping or aligning to the reference in a manner incompatible with the prescribed library type.",
            group: "Salmon Options",
        }
        range_factorization_bins: {
            description: "The range-factorization feature allows using a data-driven likelihood factorization.",
            help: "This can improve quantification accuracy on certain classes of difficult transcripts.",
            group: "Salmon Options",
        }
        fld_mean: {
            description: "Allows the user to set the expected mean fragment length of the sequencing library.",
            help: "Since the empirical fragment length distribution cannot be estimated from the mappings of single-end reads, this is only important when running Salmon with single-end reads.",
            group: "Salmon Options",
        }
        fld_sd: {
            description: "Allows the user to set the expected standard deviation of the fragment length distribution.",
            help: "Since the empirical fragment length distribution cannot be estimated from the mappings of single-end reads, this is only important when running Salmon with single-end reads.",
            group: "Salmon Options",
        }
        seq_bias: {
            description: "Passing this flag will enable it to learn and correct for sequence-specific biases in the input data.",
            group: "Salmon Options",
        }
        gc_bias: {
            description: "Passing this flag will enable it to learn and correct for fragment-level GC biases in the input data.",
            group: "Salmon Options",
        }
        pos_bias: {
            description: "Passing this flag will enable modeling of a position-specific fragment start distribution.",
            group: "Salmon Options",
        }
        use_em: {
            description: "Use the \"standard\" EM algorithm to optimize abundance estimates instead of the variational Bayesian EM algorithm.",
            group: "Salmon Options",
        }
        recover_orphans: {
            description: "This flag (which should only be used in conjunction with selective alignment), performs orphan \"rescue\" for reads.",
            group: "Salmon Options",
        }
        hard_filter: {
            description: "This flag (which should only be used with selective alignment) turns off soft filtering and range-factorized equivalence classes.",
            help: "Removes all but the equally highest scoring mappings from the equivalence class label for each fragment.",
            group: "Salmon Options",
        }
        allow_dovetail: {
            description: "Dovetailing mappings and alignments are considered discordant and discarded by default.",
            help: "If you wish to consider dovetailing mappings as concordant, you can do so by passing this flag.",
            group: "Salmon Options",
        }
        dump_eq: {
            description: "If passed, Salmon will write a file in the auxiliary directory, called eq_classes.txt.",
            help: "Contains the equivalence classes and corresponding counts that were computed during quasi-mapping.",
            group: "Salmon Options",
        }
        write_unmapped_names: {
            description: "Passing this flag will tell Salmon to write out the names of reads (or mates in paired-end reads) that do not map to the transcriptome.",
            group: "Salmon Options",
        }
        use_all_cores: {
            description: "Use all cores? Recommended for cloud environments.",
            group: "Resources",
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
        File salmon_index_tar_gz
        Array[File]+ read_one_fastqs_gz
        Array[File]? read_two_fastqs_gz
        String lib_type = "A"
        String prefix = basename(read_one_fastqs_gz[0], ".fastq.gz")
        Boolean validate_mappings = true
        Int num_bootstraps = 0
        Float incompat_prior = 0.0
        Int range_factorization_bins = 4
        Int fld_mean = 250
        Int fld_sd = 25
        Boolean seq_bias = false
        Boolean gc_bias = false
        Boolean pos_bias = false
        Boolean use_em = false
        Boolean recover_orphans = false
        Boolean hard_filter = false
        Boolean allow_dovetail = false
        Boolean dump_eq = false
        Boolean write_unmapped_names = false
        Boolean use_all_cores = false
        Int ncpu = 4
        Int modify_disk_size_gb = 0
    }

    Array[File] read_twos = select_first([read_two_fastqs_gz, []])

    Float read_one_size = size(read_one_fastqs_gz, "GB")
    Float read_two_size = size(read_twos, "GB")
    Float index_size = size(salmon_index_tar_gz, "GB")
    Int disk_size_gb = ceil((read_one_size + read_two_size + index_size) * 3) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi

        mkdir salmon_index
        tar -xzf "~{salmon_index_tar_gz}" -C salmon_index --strip-components 1

        # shellcheck disable=SC2086
        # shellcheck disable=SC2086
        salmon quant \
            -i salmon_index \
            -l "~{lib_type}" \
            ~{if length(read_twos) > 0 then "-1 " + sep(" ", squote(read_one_fastqs_gz)) + " -2 " + sep(" ", squote(read_twos)) else "-r " + sep(" ", squote(read_one_fastqs_gz))} \
            ~{if validate_mappings then "--validateMappings" else ""} \
            -p "$n_cores" \
            --numBootstraps ~{num_bootstraps} \
            --incompatPrior ~{incompat_prior} \
            --rangeFactorizationBins ~{range_factorization_bins} \
            ~{if length(read_twos) == 0 then "--fldMean " + fld_mean else ""} \
            ~{if length(read_twos) == 0 then "--fldSD " + fld_sd else ""} \
            ~{if seq_bias then "--seqBias" else ""} \
            ~{if gc_bias then "--gcBias" else ""} \
            ~{if pos_bias then "--posBias" else ""} \
            ~{if use_em then "--useEM" else ""} \
            ~{if recover_orphans then "--recoverOrphans" else ""} \
            ~{if hard_filter then "--hardFilter" else ""} \
            ~{if allow_dovetail then "--allowDovetail" else ""} \
            ~{if dump_eq then "--dumpEq" else ""} \
            ~{if write_unmapped_names then "--writeUnmappedNames" else ""} \
            -o "~{prefix}"

        tar -czf "~{prefix}.tar.gz" "~{prefix}"
    >>>

    output {
        File quant_results_tar_gz = prefix + ".tar.gz"
    }

    runtime {
        cpu: ncpu
        memory: "16 GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/salmon:1.9.0--h7e5ed60_0"
        maxRetries: 1
    }
}
