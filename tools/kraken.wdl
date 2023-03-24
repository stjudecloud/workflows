## # Kraken2
##
## Methods for bootstrapping and running [Kraken2](https://github.com/DerrickWood/kraken2)

version 1.0

task build_db {
    input {
        String db_name = "kraken2_db"
        Int kmer_len = 35
        Int minimizer_len = 31
        Int minimizer_spaces = 7
        Float load_factor = 0.7
        Int memory_gb = 8
        Int disk_size_gb = 80
        Int ncpu = 8
        Boolean detect_nproc = false
        Int max_retries = 1
    }

    String parsed_detect_nproc = if detect_nproc then "true" else ""

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if [ -n ~{parsed_detect_nproc} ]
        then
            n_cores=$(nproc)
        fi

        kraken2-build --download-taxonomy --threads "$n_cores" --db ~{db_name}

        kraken2-build --download-library archaea --threads "$n_cores" --db ~{db_name}
        kraken2-build --download-library bacteria --threads "$n_cores" --db ~{db_name}
        kraken2-build --download-library plasmid --threads "$n_cores" --db ~{db_name}
        kraken2-build --download-library viral --threads "$n_cores" --db ~{db_name}
        kraken2-build --download-library human --threads "$n_cores" --db ~{db_name}
        kraken2-build --download-library fungi --threads "$n_cores" --db ~{db_name}
        # kraken2-build --download-library plant --threads "$n_cores" --db ~{db_name}
        kraken2-build --download-library protozoa --threads "$n_cores" --db ~{db_name}
        # kraken2-build --download-library nr --threads "$n_cores" --db ~{db_name}
        # kraken2-build --download-library nt --threads "$n_cores" --db ~{db_name}
        # kraken2-build --download-library UniVec --threads "$n_cores" --db ~{db_name}
        kraken2-build --download-library UniVec_Core --threads "$n_cores" --db ~{db_name}

        kraken2-build --build \
            --kmer-len ~{kmer_len} \
            --minimizer-len ~{minimizer_len} \
            --minimizer-spaces ~{minimizer_spaces} \
            --load-factor ~{load_factor} \
            --threads "$n_cores" \
            --db ~{db_name}

        kraken2-build --clean --threads "$n_cores" --db ~{db_name}

        tar -czf "~{db_name}.tar.gz" ~{db_name}/*
    >>>
 
    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        cpu: ncpu
        docker: 'quay.io/biocontainers/kraken2:2.1.2--pl5321h9f5acd7_2'
        maxRetries: max_retries
    }

    output {
        File db = db_name + ".tar.gz"
    }

    meta {
        author: "Andrew Frantz"
        email: "andrew.frantz@stjude.org"
        description: "This WDL task builds a custom Kraken2 database."
    }
}

task kraken {
    input {
        File read1
        File read2
        File db
        String? sample_name
        Boolean store_sequences = false
        Boolean use_names = true
        Int min_base_quality = 30
        Int? memory_gb
        Int ncpu = 1
        Boolean detect_nproc = false
        Int max_retries = 1
    }

    parameter_meta {
        db: "Database for Kraken2"
    }

    String parsed_detect_nproc = if detect_nproc then "true" else ""

    Float db_size = size(db, "GiB")
    Float read1_size = size(read1, "GiB")
    Float read2_size = size(read2, "GiB")
    Int disk_size = ceil((db_size * 2) + read1_size + read2_size + 5)

    Int ram_gb = select_first([memory_gb, ceil(db_size * 1.5)])

    String inferred_basename = basename(read1, "_R1.fastq.gz")
    String sample_basename = select_first([sample_name, inferred_basename])
    String out_report = sample_basename + ".kraken2.txt"
    String out_sequences = sample_basename + ".kraken2.sequences.txt"

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if [ -n ~{parsed_detect_nproc} ]
        then
            n_cores=$(nproc)
        fi

        mkdir -p /tmp/kraken2_db/
        tar -xzf ~{db} -C /tmp/kraken2_db/ --no-same-owner

        kraken2 --db /tmp/kraken2_db/ \
            --paired \
            --output ~{if store_sequences then out_sequences else "-"} \
            --threads "$n_cores" \
            --minimum-base-quality ~{min_base_quality} \
            --report ~{out_report} \
            --report-zero-counts \
            ~{if use_names then "--use-names" else ""} \
            ~{read1} \
            ~{read2}

        rm -r /tmp/kraken2_db/
    >>>
 
    runtime {
        memory: ram_gb + " GB"
        disk: disk_size + " GB"
        cpu: ncpu
        docker: 'quay.io/biocontainers/kraken2:2.1.2--pl5321h9f5acd7_2'
        maxRetries: max_retries
    }

    output {
        File report = out_report
        File? sequences = out_sequences
    }

    meta {
        author: "Andrew Frantz"
        email: "andrew.frantz@stjude.org"
        description: "This WDL tool runs Kraken2 on a pair of fastq files."
    }
}
