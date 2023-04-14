## # Kraken2
##
## Methods for bootstrapping and running [Kraken2](https://github.com/DerrickWood/kraken2)

version 1.0

task download_taxonomy {
    input {
        Int memory_gb = 4
        Int disk_size_gb = 60
        Int max_retries = 3
    }

    parameter_meta {
        memory_gb: "RAM to allocate for task"
        disk_size_gb: "Disk space to allocate for task"
        max_retries: "Number of times to retry in case of failure"
    }

    String db_name = "kraken2_taxonomy"

    command <<<
        set -euo pipefail

        kraken2-build --download-taxonomy \
            --use-ftp \
            --db ~{db_name} 2>&1 \
            | awk '/gunzip:/ { print; exit 42 } !/gunzip:/ { print }' 1>&2

        tar -C ~{db_name}/ -czf "~{db_name}.tar.gz" .

        rm -r ~{db_name}
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
        description: "This WDL task downloads the NCBI taxonomy which Kraken2 uses to create a tree and taxon map during the database build."
    }
}

task download_library {
    input {
        String library_name
        Int memory_gb = 4
        Int added_disk_size_gb = 0
        Int max_retries = 3
    }

    String db_name = "kraken2_"+library_name+"_library"

    parameter_meta {
        library_name: {
            description: "Library to download. Note that due to the large size of the `nt` library, the `added_disk_size_gb` parameter might have to be used, depending on your backend. This option was not tested, but we estimate the size required to be around 2TB."
            choices: ['archaea', 'bacteria', 'plasmid', 'viral', 'human', 'fungi', 'plant', 'protozoa', 'nt', 'UniVec', 'UniVec_Core']
        }
        memory_gb: "RAM to allocate for task"
        added_disk_size_gb: "Additional disk space to allocate for task. Default disk size is determined dynamically based on `library_name`."
        max_retries: "Number of times to retry in case of failure"
    }

    Int disk_size_gb = (if library_name=="bacteria" then 300 else 20) + added_disk_size_gb

    command <<<
        set -euo pipefail

        kraken2-build --download-library \
            ~{library_name} \
            --use-ftp \
            --db ~{db_name} 2>&1 \
            | awk '/gunzip:/ { print; exit 42 } !/gunzip:/ { print }' 1>&2

        tar -C ~{db_name}/ -czf "~{db_name}.tar.gz" .

        rm -r ~{db_name}
    >>>
 
    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        cpu: 1
        docker: 'ghcr.io/stjudecloud/kraken2:branch-kraken2-1.0.0'
        maxRetries: max_retries
    }

    output {
        File library = db_name + ".tar.gz"
    }

    meta {
        author: "Andrew Frantz"
        email: "andrew.frantz@stjude.org"
        description: "This WDL task downloads a predefined library of reference genomes from NCBI. Detailed organism list for libraries (except nt) available at: https://ftp.ncbi.nlm.nih.gov/genomes/refseq/."
    }
}

task create_library_from_fastas {
    input {
        Array[File] fastas
        Int memory_gb = 4
        Int added_disk_size_gb = 0
        Int max_retries = 1
    }

    String db_name = "kraken2_custom_library"

    parameter_meta {
        fastas: "Array of gzipped FASTA files. Each FASTA sequence ID must contain either an NCBI accession number or an explicit assignment of the taxonomy ID using `kraken:taxid`"
        memory_gb: "RAM to allocate for task"
        added_disk_size_gb: "Additional disk space to allocate for task. Default disk size is determined dynamically based on `fastas` size."
        max_retries: "Number of times to retry in case of failure"
    }

    Int fastas_size = ceil(size(fastas, "GiB"))
    Int disk_size_gb = fastas_size * 5 + added_disk_size_gb

    command <<<
        set -euo pipefail

        >&2 echo "*** start adding custom FASTAs ***"
        echo "~{sep="\n" fastas}" > fastas.txt
        while read -r fasta; do
            gunzip -c "$fasta" > tmp.fa
            kraken2-build --add-to-library tmp.fa --db ~{db_name}
        done < fastas.txt
        rm tmp.fa
        >&2 echo "*** done adding custom FASTAs ***"

        tar -C ~{db_name}/ -czf "~{db_name}.tar.gz" .

        rm -r ~{db_name}
    >>>
 
    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        cpu: 1
        docker: 'quay.io/biocontainers/kraken2:2.1.2--pl5321h9f5acd7_2'
        maxRetries: max_retries
    }

    output {
        File custom_library = db_name + ".tar.gz"
    }

    meta {
        author: "Andrew Frantz"
        email: "andrew.frantz@stjude.org"
        description: "This WDL task adds custom entries from FASTA files to a Kraken2 DB."
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
        Int added_memory_gb = 0
        Int added_disk_size_gb = 0
        Int ncpu = 1
        Boolean detect_nproc = false
        Int max_retries = 1
    }

    parameter_meta {
        tarballs: "Tarballs containing the NCBI taxonomy (generated by the `download_taxonomy` task) and at least one library (generated by the `download_library` task). Tarballs must not have a root directory."
        db_name: "Name for output in compressed, archived format. The suffix `.tar.gz` will be added."
        kmer_len: "K-mer length in bp that will be used to build the database"
        minimizer_len: "Minimizer length in bp that will be used to build the database"
        minimizer_spaces: "Number of characters in minimizer that are ignored in comparisons"
        max_db_size_gb: "Maximum number of GBs for Kraken 2 hash table; if the Kraken 2 estimator determines more would normally be needed, the reference library will be downsampled to fit."
        added_memory_gb: "Additional RAM to allocate for task. Default RAM is allocated dynamically based on the database size."
        added_disk_size_gb: "Additional disk space to allocate for task. Default disk size is determined dynamically based on size of the input `tarballs`."
        ncpu: "Number of cores to allocate for task"
        detect_nproc: "Use all available cores"
        max_retries: "Number of times to retry in case of failure"
    }

    Int tarballs_size = ceil(size(tarballs, "GiB"))
    Int disk_size_gb = tarballs_size * 6 + added_disk_size_gb
    Int memory_gb = (if (max_db_size_gb > 0) then (ceil(max_db_size_gb * 1.2)) else (tarballs_size * 2)) + added_memory_gb

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if [ "~{detect_nproc}" = "true" ]; then
            n_cores=$(nproc)
        fi

        >&2 echo "*** start unpacking tarballs ***"
        echo "~{sep="\n" tarballs}" > tarballs.txt
        mkdir ~{db_name}
        while read -r tarball; do
            tar -xzf "$tarball" -C ~{db_name} --no-same-owner
        done < tarballs.txt
        >&2 echo "*** done unpacking tarballs ***"

        >&2 echo "*** start DB build ***"
        kraken2-build --build \
            --kmer-len ~{kmer_len} \
            --minimizer-len ~{minimizer_len} \
            --minimizer-spaces ~{minimizer_spaces} \
            ~{if (max_db_size_gb > 0) then "--max-db-size" else ""} ~{if (max_db_size_gb > 0) then max_db_size_gb + "000000000" else ""} \
            --threads "$n_cores" \
            --db ~{db_name}

        >&2 echo "*** start DB clean ***"
        kraken2-build --clean --threads "$n_cores" --db ~{db_name}
        >&2 echo "*** done ***"

        >&2 echo "*** tarballing DB ***"
        tar -C ~{db_name}/ -czf "~{db_name}.tar.gz" .

        rm -r ~{db_name}
    >>>
 
    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        cpu: ncpu
        docker: 'quay.io/biocontainers/kraken2:2.1.2--pl5321h9f5acd7_2'
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
        read1: "Gzipped FastQ file with 1st reads in pair"
        read2: "Gzipped FastQ file with 2nd reads in pair"
        db: "Kraken2 database. Can be generated with `make-qc-reference.wdl`. Must be a flat tarball without a root directory."
        sample_name: "Name for sample. If missing will be inferred by removing the suffix '_R1.fastq.gz' from the `read1` filename."
        store_sequences: "Store and output main Kraken2 output in addition to the summary report"
        use_names: "Print scientific names instead of just taxids"
        min_base_quality: "Minimum base quality used in classification"
        memory_gb: "RAM to allocate for task. If missing will be dynamically allocated based on database size."
        ncpu: "Number of cores to allocate for task"
        detect_nproc: "Use all available cores"
        max_retries: "Number of times to retry in case of failure"
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
