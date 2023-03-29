## # Kraken2
##
## Methods for bootstrapping and running [Kraken2](https://github.com/DerrickWood/kraken2)

version 1.0

task build_db_full {
    input {
        String db_name = "kraken2_db"
        File? base_db
        Array[File] fastas = []
        Int kmer_len = 35
        Int minimizer_len = 31
        Int minimizer_spaces = 7
        Int max_db_size_gb = -1
        Float load_factor = 0.7
        Int memory_gb = 96
        Int added_disk_size_gb = 200
        Int ncpu = 1
        Boolean detect_nproc = false
        Int max_retries = 1
    }

    Int fastas_size = ceil(size(fastas, "GiB"))
    Int disk_size_gb = fastas_size * 2 + added_disk_size_gb

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if [ "~{detect_nproc}" = "true" ]; then
            n_cores=$(nproc)
        fi

        if [ ~{defined(base_db)} = "true" ]; then
            >&2 echo "*** start unpacking base DB ***"
            mkdir ~{db_name}
            tar -xzf ~{base_db} -C ~{db_name} --no-same-owner
            >&2 echo "*** done unpacking base DB ***"
        fi

        >&2 echo "*** start downloading taxonomy ***"
        kraken2-build --download-taxonomy --use-ftp --threads "$n_cores" --db ~{db_name}
        >&2 echo "*** done downloading taxonomy ***"

        if [ "~{length(fastas) > 0}" = "true" ]; then
            >&2 echo "*** start adding custom FASTAs ***"
            echo "~{sep="\n" fastas}" > fastas.txt
            while read -r fasta; do
                gunzip -c "$fasta" > tmp.fa
                kraken2-build --add-to-library tmp.fa --use-ftp --threads "$n_cores" --db ~{db_name}
            done < fastas.txt
            rm tmp.fa
            >&2 echo "*** done adding custom FASTAs ***"
        fi

        if [ ~{defined(base_db)} = "false" ]; then
            >&2 echo "*** start downloading archaea ***"
            kraken2-build --download-library archaea --use-ftp --threads "$n_cores" --db ~{db_name}
            >&2 echo "*** start downloading bacteria ***"
            kraken2-build --download-library bacteria --use-ftp --threads "$n_cores" --db ~{db_name}
            >&2 echo "*** start downloading plasmid ***"
            kraken2-build --download-library plasmid --use-ftp --threads "$n_cores" --db ~{db_name}
            >&2 echo "*** start downloading viral ***"
            kraken2-build --download-library viral --use-ftp --threads "$n_cores" --db ~{db_name}
            >&2 echo "*** start downloading human ***"
            kraken2-build --download-library human --use-ftp --threads "$n_cores" --db ~{db_name}
            >&2 echo "*** start downloading fungi ***"
            kraken2-build --download-library fungi --use-ftp --threads "$n_cores" --db ~{db_name}
            # kraken2-build --download-library plant --use-ftp --threads "$n_cores" --db ~{db_name}
            >&2 echo "*** start downloading protozoa ***"
            kraken2-build --download-library protozoa --use-ftp --threads "$n_cores" --db ~{db_name}
            # kraken2-build --download-library nr --use-ftp --threads "$n_cores" --db ~{db_name}
            # kraken2-build --download-library nt --use-ftp --threads "$n_cores" --db ~{db_name}
            # kraken2-build --download-library UniVec --use-ftp --threads "$n_cores" --db ~{db_name}
            >&2 echo "*** start downloading UniVec_Core ***"
            kraken2-build --download-library UniVec_Core --use-ftp --threads "$n_cores" --db ~{db_name}
        fi

        >&2 echo "*** start DB build ***"
        kraken2-build --build \
            --kmer-len ~{kmer_len} \
            --minimizer-len ~{minimizer_len} \
            --minimizer-spaces ~{minimizer_spaces} \
            ~{if (max_db_size_gb > 0) then "--max-db-size" else ""} ~{if (max_db_size_gb > 0) then max_db_size_gb + "000000000" else ""} \
            --load-factor ~{load_factor} \
            --threads "$n_cores" \
            --db ~{db_name}

        >&2 echo "*** start DB clean ***"
        kraken2-build --clean --threads "$n_cores" --db ~{db_name}
        >&2 echo "*** done ***"

        >&2 echo "*** tarballing DB ***"
        tar -czf "~{db_name}.tar.gz" ~{db_name}/*
    >>>
 
    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        cpu: ncpu
        docker: 'ghcr.io/stjudecloud/kraken2:branch-kraken2-1.0.0'
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

task download_taxonomy {
    input {
        String db_name = "kraken2_db"
        Int memory_gb = 4
        Int disk_size_gb = 60
        Int max_retries = 1
    }

    command <<<
        set -euo pipefail

        kraken2-build --download-taxonomy --use-ftp --db ~{db_name}

        tar -czf "~{db_name}.tar.gz" ~{db_name}/*
    >>>
 
    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        cpu: 1
        docker: 'ghcr.io/stjudecloud/kraken2:branch-kraken2-1.0.0'
        maxRetries: max_retries
    }

    output {
        File taxonomy = db_name + ".tar.gz"
    }

    meta {
        author: "Andrew Frantz"
        email: "andrew.frantz@stjude.org"
        description: "This WDL task downloads the Kraken2 taxonomy."
    }
}

task download_library {
    input {
        String library
        String db_name = "kraken2_db"
        Int memory_gb = 4
        Int disk_size_gb = 60
        Int max_retries = 1
    }

    command <<<
        set -euo pipefail

        kraken2-build --download-library ~{library} --use-ftp --db ~{db_name}

        tar -czf "~{db_name}.tar.gz" ~{db_name}/*
    >>>
 
    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        cpu: 1
        docker: 'ghcr.io/stjudecloud/kraken2:branch-kraken2-1.0.0'
        maxRetries: max_retries
    }

    output {
        File db = db_name + ".tar.gz"
    }

    meta {
        author: "Andrew Frantz"
        email: "andrew.frantz@stjude.org"
        description: "This WDL task downloads a Kraken2 library."
    }
}

task add_custom_fastas_to_db {
    input {
        Array[File] fastas
        String db_name = "kraken2_db"
        Int memory_gb = 4
        Int disk_size_gb = 60
        Int max_retries = 1
    }

    command <<<
        set -euo pipefail

        >&2 echo "*** start adding custom FASTAs ***"
        echo "~{sep="\n" fastas}" > fastas.txt
        while read -r fasta; do
            gunzip -c "$fasta" > tmp.fa
            kraken2-build --add-to-library tmp.fa --use-ftp --db ~{db_name}
        done < fastas.txt
        rm tmp.fa
        >&2 echo "*** done adding custom FASTAs ***"

        tar -czf "~{db_name}.tar.gz" ~{db_name}/*
    >>>
 
    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        cpu: 1
        docker: 'ghcr.io/stjudecloud/kraken2:branch-kraken2-1.0.0'
        maxRetries: max_retries
    }

    output {
        File db = db_name + ".tar.gz"
    }

    meta {
        author: "Andrew Frantz"
        email: "andrew.frantz@stjude.org"
        description: "This WDL task adds custom FASTAs to a Kraken2 DB."
    }
}

task build_db {
    input {
        Array[File] tarballs
        String db_name = "kraken2_db"
        Int kmer_len = 35
        Int minimizer_len = 31
        Int minimizer_spaces = 7
        Int max_db_size_gb = -1
        Float load_factor = 0.7
        Int memory_gb = 96
        Int disk_size_gb = 200
        Int ncpu = 1
        Boolean detect_nproc = false
        Int max_retries = 1
    }

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if [ "~{detect_nproc}" = "true" ]; then
            n_cores=$(nproc)
        fi

        >&2 echo "*** start unpacking tarballs ***"
        echo "~{sep="\n" tarballs}" > tarballs.txt
        while read -r tarball; do
            mkdir ~{db_name}
            tar -xzf "$tarball" -C ~{db_name} --no-same-owner
        done < tarballs.txt
        >&2 echo "*** done unpacking tarballs ***"

        >&2 echo "*** start DB build ***"
        kraken2-build --build \
            --kmer-len ~{kmer_len} \
            --minimizer-len ~{minimizer_len} \
            --minimizer-spaces ~{minimizer_spaces} \
            ~{if (max_db_size_gb > 0) then "--max-db-size" else ""} ~{if (max_db_size_gb > 0) then max_db_size_gb + "000000000" else ""} \
            --load-factor ~{load_factor} \
            --threads "$n_cores" \
            --db ~{db_name}

        >&2 echo "*** start DB clean ***"
        kraken2-build --clean --threads "$n_cores" --db ~{db_name}
        >&2 echo "*** done ***"

        >&2 echo "*** tarballing DB ***"
        tar -czf "~{db_name}.tar.gz" ~{db_name}/*
    >>>
 
    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        cpu: ncpu
        docker: 'ghcr.io/stjudecloud/kraken2:branch-kraken2-1.0.0'
        maxRetries: max_retries
    }

    output {
        File built_db = db_name + ".tar.gz"
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
        Int min_base_quality = 0
        Int? memory_gb
        Int ncpu = 1
        Boolean detect_nproc = false
        Int max_retries = 1
    }

    parameter_meta {
        db: "Database for Kraken2"
    }

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
        if [ "~{detect_nproc}" = "true" ]; then
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
