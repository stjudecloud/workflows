## # Kraken2
##
## Methods for bootstrapping and running [Kraken2](https://github.com/DerrickWood/kraken2)

version 1.0

task download_taxonomy {
    input {
        Boolean protein = false
        Int memory_gb = 4
        Int disk_size_gb = 60
        Int max_retries = 3
    }

    parameter_meta {
        protein: "Construct a protein database?"
        memory_gb: "RAM to allocate for task"
        disk_size_gb: "Disk space to allocate for task"
        max_retries: "Number of times to retry in case of failure"
    }

    String db_name = "kraken2_taxonomy"

    command <<<
        set -euo pipefail

        kraken2-build --download-taxonomy \
            ~{if protein then "--protein" else ""} \
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
        docker: 'ghcr.io/stjudecloud/kraken2:2.1.2-0'
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
        Boolean protein = false
        Int memory_gb = 4
        Int added_disk_size_gb = 0
        Int max_retries = 3
    }

    String db_name = "kraken2_"+library_name+"_library"

    parameter_meta {
        library_name: {
            description: "Library to download. Note that `protein` must equal `true` if downloading the `nr` library, and `protein` must equal `false` if downloading the `UniVec` or `UniVec_Core` library."
            choices: ['archaea', 'bacteria', 'plasmid', 'viral', 'human', 'fungi', 'plant', 'protozoa', 'nt', 'nr', 'UniVec', 'UniVec_Core']
        }
        protein: "Construct a protein database?"
        memory_gb: "RAM to allocate for task"
        added_disk_size_gb: "Additional disk space to allocate for task. Default disk size is determined dynamically based on `library_name`. Note that the default sizes are adequate as of April 2023, but new genomes are constantly being added to the NCBI database. More disk space may be required depending on when in the future this task is run."
        max_retries: "Number of times to retry in case of failure"
    }

    Int disk_size_gb = (
        if library_name=="bacteria" then 300
        else if library_name=="nr" then 600
        else if library_name=="nt" then 2500
        else 20
    ) + added_disk_size_gb

    command <<<
        set -euo pipefail

        kraken2-build --download-library \
            ~{library_name} \
            ~{if protein then "--protein" else ""} \
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
        docker: 'ghcr.io/stjudecloud/kraken2:2.1.2-0'
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
        Boolean protein = false
        Int memory_gb = 4
        Int added_disk_size_gb = 0
        Int max_retries = 1
    }

    String db_name = "kraken2_custom_library"

    parameter_meta {
        fastas: "Array of gzipped FASTA files. Each FASTA sequence ID must contain either an NCBI accession number or an explicit assignment of the taxonomy ID using `kraken:taxid`"
        protein: "Construct a protein database?"
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
            kraken2-build \
                ~{if protein then "--protein" else ""} \
                --add-to-library tmp.fa \
                --db ~{db_name}
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
        Boolean protein = false
        Int kmer_len = if protein then 15 else 35
        Int minimizer_len = if protein then 12 else 31
        Int minimizer_spaces = if protein then 0 else 7
        Int max_db_size_gb = -1
        Int added_memory_gb = 0
        Int added_disk_size_gb = 0
        Int ncpu = 1
        Boolean use_all_cores = false
        Int max_retries = 1
    }

    parameter_meta {
        tarballs: "Tarballs containing the NCBI taxonomy (generated by the `download_taxonomy` task) and at least one library (generated by the `download_library` task). Tarballs must not have a root directory."
        db_name: "Name for output in compressed, archived format. The suffix `.tar.gz` will be added."
        protein: "Construct a protein database?"
        kmer_len: "K-mer length in bp that will be used to build the database"
        minimizer_len: "Minimizer length in bp that will be used to build the database"
        minimizer_spaces: "Number of characters in minimizer that are ignored in comparisons"
        max_db_size_gb: "Maximum number of GBs for Kraken 2 hash table; if the Kraken 2 estimator determines more would normally be needed, the reference library will be downsampled to fit."
        added_memory_gb: "Additional RAM to allocate for task. Default RAM is allocated dynamically based on the database size."
        added_disk_size_gb: "Additional disk space to allocate for task. Default disk size is determined dynamically based on size of the input `tarballs`."
        ncpu: "Number of cores to allocate for task"
        use_all_cores: "Use all cores. Recommended for cloud environments. Not recommended for cluster environments."
        max_retries: "Number of times to retry in case of failure"
    }

    Int tarballs_size = ceil(size(tarballs, "GiB"))
    Int disk_size_gb = tarballs_size * 6 + added_disk_size_gb
    Int memory_gb = (
        if (max_db_size_gb > 0)
        then (ceil(max_db_size_gb * 1.2))
        else (tarballs_size * 2)
    ) + added_memory_gb

    String max_db_size_bytes = max_db_size_gb + "000000000"

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if [ "~{use_all_cores}" = "true" ]; then
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
            ~{if protein then "--protein" else ""} \
            --kmer-len ~{kmer_len} \
            --minimizer-len ~{minimizer_len} \
            --minimizer-spaces ~{minimizer_spaces} \
            ~{if (max_db_size_gb > 0) then "--max-db-size" else ""} ~{if (max_db_size_gb > 0) then max_db_size_bytes else ""} \
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
        String? prefix = basename(read1, "_R1.fastq.gz")
        Boolean store_sequences = false
        Boolean use_names = true
        Int min_base_quality = 0
        Int modify_memory_gb = 0
        Int modify_disk_size_gb = 0
        Int ncpu = 1
        Boolean use_all_cores = false
        Int max_retries = 1
    }

    parameter_meta {
        read1: "Gzipped FastQ file with 1st reads in pair"
        read2: "Gzipped FastQ file with 2nd reads in pair"
        db: "Kraken2 database. Can be generated with `make-qc-reference.wdl`. Must be a tarball without a root directory."
        sample_name: "Name for sample. If missing will be inferred by removing the suffix '_R1.fastq.gz' from the `read1` filename."
        store_sequences: "Store and output main Kraken2 output in addition to the summary report"
        use_names: "Print scientific names instead of just taxids"
        min_base_quality: "Minimum base quality used in classification"
        memory_gb: "RAM to allocate for task. If missing will be dynamically allocated based on database size."
        ncpu: "Number of cores to allocate for task"
        use_all_cores: "Use all cores? Recommended for cloud environments. Not recommended for cluster environments."
        max_retries: "Number of times to retry in case of failure"
    }

    Float db_size = size(db, "GiB")
    Float read1_size = size(read1, "GiB")
    Float read2_size = size(read2, "GiB")
    Int disk_size_gb_calculation = ceil(
        (db_size * 2) + read1_size + read2_size
    ) + modify_disk_size_gb
    Int disk_size_gb = if store_sequences
        then disk_size_gb_calculation + ceil(read1_size + read2_size)
        else disk_size_gb_calculation

    Int memory_gb = ceil(db_size * 2) + modify_memory_gb

    String out_report = prefix + ".kraken2.txt"
    String out_sequences = prefix + ".kraken2.sequences.txt"

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if [ "~{use_all_cores}" = "true" ]; then
            n_cores=$(nproc)
        fi

        mkdir -p kraken2_db/
        tar -xzf ~{db} -C kraken2_db/ --no-same-owner

        kraken2 --db kraken2_db/ \
            --paired \
            --output ~{if store_sequences then out_sequences else "-"} \
            --threads "$n_cores" \
            --minimum-base-quality ~{min_base_quality} \
            --report ~{out_report} \
            --report-zero-counts \
            ~{if use_names then "--use-names" else ""} \
            ~{read1} \
            ~{read2}

        if [ "~{store_sequences}" = "true" ]; then
            gzip ~{out_sequences}
        fi

        rm -r kraken2_db/
    >>>
 
    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        cpu: ncpu
        docker: 'quay.io/biocontainers/kraken2:2.1.2--pl5321h9f5acd7_2'
        maxRetries: max_retries
    }

    output {
        File report = out_report
        File? sequences = out_sequences + ".gz"
    }

    meta {
        author: "Andrew Frantz"
        email: "andrew.frantz@stjude.org"
        description: "This WDL tool runs Kraken2 on a pair of fastq files."
    }
}
