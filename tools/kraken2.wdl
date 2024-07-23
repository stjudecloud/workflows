## [Homepage](https://github.com/DerrickWood/kraken2)

version 1.1

task download_taxonomy {
    meta {
        description: "Downloads the NCBI taxonomy which Kraken2 uses to create a tree and taxon map during the database build"
        outputs: {
            taxonomy: "The NCBI taxonomy, which is needed by the `build_db` task. This output is not human-readable or meant for anything other than building a Kraken2 database."
        }
    }

    parameter_meta {
        protein: "Construct a protein database?"
    }

    input {
        Boolean protein = false
    }

    String db_name = "kraken2_taxonomy"

    command <<<
        set -euo pipefail

        kraken2-build --download-taxonomy \
            ~{if protein then "--protein" else ""} \
            --use-ftp \
            --db ~{db_name} 2>&1 \
            | awk '/gunzip:/ { print; exit 42 } !/gunzip:/ { print }' 1>&2

        tar -C ~{db_name} -czf "~{db_name}.tar.gz" .

        rm -r ~{db_name}
    >>>

    output {
        File taxonomy = db_name + ".tar.gz"
    }

    runtime {
        memory: "4 GB"
        disks: "60 GB"
        container: "quay.io/biocontainers/kraken2:2.1.3--pl5321hdcf5f25_0"
        maxRetries: 3
    }
}

task download_library {
    meta {
        description: "Downloads a predefined library of reference genomes from NCBI. Detailed organism list for libraries (except nt) available [here](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/)"
        outputs: {
            library: "A library of reference genomes, which is needed by the `build_db` task. This output is not human-readable or meant for anything other than building a Kraken2 database."
        }
    }

    parameter_meta {
        library_name: {
            description: "Library to download. Note that `protein` must equal `true` if downloading the `nr` library, and `protein` must equal `false` if downloading the `UniVec` or `UniVec_Core` library.",
            choices: [
                "archaea",
                "bacteria",
                "plasmid",
                "viral",
                "human",
                "fungi",
                "plant",
                "protozoa",
                "nt",
                "nr",
                "UniVec",
                "UniVec_Core",
            ],
        }
        protein: "Construct a protein database?"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation, specified in GB. Default disk size is determined dynamically based on `library_name`. Note that the default sizes are adequate as of April 2023, but new genomes are constantly being added to the NCBI database. More disk space may be required depending on when in the future this task is run."
    }

    input {
        String library_name
        Boolean protein = false
        Int modify_disk_size_gb = 0
    }

    String db_name = "kraken2_"+library_name+"_library"

    Int disk_size_gb = (
        (
            if library_name=="bacteria" then 300
            else if library_name=="nr" then 600
            else if library_name=="nt" then 2500
            else 25
        ) + modify_disk_size_gb
    )

    command <<<
        set -euo pipefail

        kraken2-build --download-library \
            ~{library_name} \
            ~{if protein then "--protein" else ""} \
            --use-ftp \
            --db ~{db_name} 2>&1 \
            | awk '/gunzip:/ { print; exit 42 } !/gunzip:/ { print }' 1>&2

        tar -C ~{db_name} -czf "~{db_name}.tar.gz" .

        rm -r ~{db_name}
    >>>

    output {
        File library = db_name + ".tar.gz"
    }

    runtime {
        memory: "4 GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/kraken2:2.1.3--pl5321hdcf5f25_0"
        maxRetries: 3
    }
}

task create_library_from_fastas {
    meta {
        description: "Adds custom entries from FASTA files to a Kraken2 DB"
        outputs: {
            custom_library: "Kraken2 compatible library, which is needed by the `build_db` task. This output is not human-readable or meant for anything other than building a Kraken2 database."
        }
    }

    parameter_meta {
        fastas_gz: "Array of gzipped FASTA files. Each FASTA sequence ID must contain either an NCBI accession number or an explicit assignment of the taxonomy ID using `kraken:taxid`"
        protein: "Construct a protein database?"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        Array[File] fastas_gz
        Boolean protein = false
        Int modify_disk_size_gb = 0
    }

    String db_name = "kraken2_custom_library"

    Float fastas_size = size(fastas_gz, "GiB")
    Int disk_size_gb = ceil(fastas_size * 5) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        >&2 echo "*** start adding custom FASTAs ***"
        echo "~{sep("\n", fastas_gz)}" > fastas.txt
        while read -r fasta; do
            gunzip -c "$fasta" > tmp.fa
            kraken2-build \
                ~{if protein then "--protein" else ""} \
                --add-to-library tmp.fa \
                --db ~{db_name}
        done < fastas.txt
        rm tmp.fa
        >&2 echo "*** done adding custom FASTAs ***"

        tar -C ~{db_name} -czf "~{db_name}.tar.gz" .

        rm -r ~{db_name}
    >>>

    output {
        File custom_library = db_name + ".tar.gz"
    }

    runtime {
        memory: "4 GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/kraken2:2.1.3--pl5321hdcf5f25_0"
        maxRetries: 1
    }
}

task build_db {
    meta {
        description: "Builds a custom Kraken2 database"
        outputs: {
            built_db: "A complete Kraken2 database"
        }
    }

    parameter_meta {
        tarballs: "Tarballs containing the NCBI taxonomy (generated by the `download_taxonomy` task) and at least one library (generated by the `download_library` or `create_library_from_fastas` task). Tarballs must not have a root directory."
        db_name: {
            description: "Name for output in compressed, archived format. The suffix `.tar.gz` will be added.",
            common: true,
        }
        protein: "Construct a protein database?"
        use_all_cores: {
            description: "Use all cores? Recommended for cloud environments.",
            common: true,
        }
        kmer_len: "K-mer length in bp that will be used to build the database"
        minimizer_len: "Minimizer length in bp that will be used to build the database"
        minimizer_spaces: "Number of characters in minimizer that are ignored in comparisons"
        max_db_size_gb: "Maximum number of GBs for Kraken 2 hash table; if the Kraken 2 estimator determines more would normally be needed, the reference library will be downsampled to fit."
        ncpu: {
            description: "Number of cores to allocate for task",
            common: true,
        }
        modify_memory_gb: "Add to or subtract from dynamic memory allocation. Default memory is determined by the size of the inputs. Specified in GB."
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        Array[File] tarballs
        String db_name = "kraken2_db"
        Boolean protein = false
        Boolean use_all_cores = false
        Int kmer_len = if protein then 15 else 35
        Int minimizer_len = if protein then 12 else 31
        Int minimizer_spaces = if protein then 0 else 7
        Int max_db_size_gb = -1
        Int ncpu = 4
        Int modify_memory_gb = 0
        Int modify_disk_size_gb = 0
    }

    Float tarballs_size = size(tarballs, "GiB")
    Int disk_size_gb = ceil(tarballs_size * 6) + 10 + modify_disk_size_gb
    Int memory_gb = (
        (
            if (max_db_size_gb > 0)
            then ceil(max_db_size_gb * 1.2)
            else ceil(tarballs_size * 2)
        ) + modify_memory_gb
    )

    String max_db_size_bytes = "~{max_db_size_gb}000000000"

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi

        >&2 echo "*** start unpacking tarballs ***"
        echo "~{sep("\n", tarballs)}" > tarballs.txt
        mkdir ~{db_name}
        while read -r tarball; do
            tar -xzf "$tarball" -C ~{db_name} --no-same-owner
        done < tarballs.txt
        >&2 echo "*** done unpacking tarballs ***"

        >&2 echo "*** start DB build ***"
        kraken2-build --build \
            ~{if protein then "--protein" else ""} \
            --kmer-len ~{kmer_len} \
            --minimizer-len ~{minimizer_len} \
            --minimizer-spaces ~{minimizer_spaces} \
            ~{if (max_db_size_gb > 0)
                then "--max-db-size " + max_db_size_bytes
                else ""
            } \
            --threads "$n_cores" \
            --db ~{db_name}

        >&2 echo "*** start DB clean ***"
        kraken2-build --clean --threads "$n_cores" --db ~{db_name}
        >&2 echo "*** done ***"

        >&2 echo "*** tarballing DB ***"
        tar -C ~{db_name}/ -czf "~{db_name}.tar.gz" .

        rm -r ~{db_name}
    >>>

    output {
        File built_db = db_name + ".tar.gz"
    }

    runtime {
        cpu: ncpu
        memory: "~{memory_gb} GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/kraken2:2.1.3--pl5321hdcf5f25_0"
        maxRetries: 1
    }
}

task kraken {
    meta {
        description: "Runs Kraken2 on a pair of fastq files"
        outputs: {
            report: {
                description: "A Kraken2 summary report",
                external_help: "https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#sample-report-output-format",
            },
            sequences: {
                description: "Detailed Kraken2 output that has been gzipped",
                external_help: "https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#standard-kraken-output-format",
            },
        }
    }

    parameter_meta {
        read_one_fastq_gz: "Gzipped FASTQ file with 1st reads in pair"
        read_two_fastq_gz: "Gzipped FASTQ file with 2nd reads in pair"
        db: "Kraken2 database. Can be generated with `make-qc-reference.wdl`. Must be a tarball without a root directory."
        prefix: "Prefix for the Kraken2 output files. The extensions `.kraken2.txt` and `.kraken2.sequences.txt.gz` will be added."
        store_sequences: {
            description: "Store and output main Kraken2 output in addition to the summary report?",
            common: true,
        }
        use_names: "Print scientific names instead of just taxids?"
        use_all_cores: {
            description: "Use all cores? Recommended for cloud environments.",
            common: true,
        }
        min_base_quality: "Minimum base quality used in classification"
        ncpu: {
            description: "Number of cores to allocate for task",
            common: true,
        }
        modify_memory_gb: "Add to or subtract from dynamic memory allocation. Default memory is determined by the size of the inputs. Specified in GB."
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File read_one_fastq_gz
        File read_two_fastq_gz
        File db
        String prefix = sub(
            basename(read_one_fastq_gz),
            "([_\\.][rR][12])?(\\.subsampled)?\\.(fastq|fq)(\\.gz)?$",
            ""
        )
        Boolean store_sequences = false
        Boolean use_names = true
        Boolean use_all_cores = false
        Int min_base_quality = 0
        Int ncpu = 4
        Int modify_memory_gb = 0
        Int modify_disk_size_gb = 0
    }

    Float db_size = size(db, "GiB")
    Float read1_size = size(read_one_fastq_gz, "GiB")
    Float read2_size = size(read_two_fastq_gz, "GiB")
    Int disk_size_gb_calculation = (
        ceil((db_size * 2) + read1_size + read2_size) + 10 + modify_disk_size_gb
    )
    Int disk_size_gb = if store_sequences
        then disk_size_gb_calculation + ceil(read1_size + read2_size)
        else disk_size_gb_calculation

    Int memory_gb = ceil(db_size * 2) + modify_memory_gb

    String out_report = prefix + ".kraken2.txt"
    String out_sequences = prefix + ".kraken2.sequences.txt"

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi

        mkdir kraken2_db/
        tar -xzf ~{db} -C kraken2_db/ --no-same-owner

        kraken2 --db kraken2_db/ \
            --paired \
            --output ~{if store_sequences then out_sequences else "-"} \
            --threads "$n_cores" \
            --minimum-base-quality ~{min_base_quality} \
            --report ~{out_report} \
            --report-zero-counts \
            ~{if use_names then "--use-names" else ""} \
            ~{read_one_fastq_gz} \
            ~{read_two_fastq_gz}

        if ~{store_sequences}; then
            gzip ~{out_sequences}
        fi

        rm -r kraken2_db/
    >>>

    output {
        File report = out_report
        File? sequences = out_sequences + ".gz"
    }

    runtime {
        cpu: ncpu
        memory: "~{memory_gb} GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/kraken2:2.1.3--pl5321hdcf5f25_0"
        maxRetries: 1
    }
}
