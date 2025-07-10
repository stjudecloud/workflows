version 1.1

import "../../tools/kraken2.wdl"
import "../../tools/util.wdl"

workflow qc_reference {
    meta {
        name: "Quality Check Reference"
        description: "Downloads and creates all reference files needed to run the `quality_check_standard` workflow"
        warning: "See `kraken2.download_library.meta.warning` for information regarding common failures."
        category: "Reference"
        outputs: {
            reference_fa: "FASTA format reference file",
            gtf: "GTF feature file",
            kraken_db: "A complete Kraken2 database",
            coverage_beds: "BED file for each feature type to use for coverage calculation",
        }
        allowNestedInputs: true
    }

    parameter_meta {
        reference_fa_url: "URL to retrieve the reference FASTA file from"
        reference_fa_name: "Name of output reference FASTA file"
        gtf_url: "URL to retrieve the reference GTF file from"
        gtf_name: "Name of output GTF file"
        kraken_fastas: "Array of gzipped FASTA files. Each sequence's ID must contain either an NCBI accession number or an explicit assignment of the taxonomy ID using `kraken:taxid`"
        kraken_fasta_urls: "URLs for any additional FASTA files in NCBI format to download and include in the Kraken2 database. This allows the addition of individual genomes (or other sequences) of interest."
        kraken_libraries: {
            description: "List of kraken libraries to download",
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
                "UniVec",
                "UniVec_Core",
            ],
        }
        coverage_feature_types: {
            description: "List of feature types to use for coverage calculation",
            help: "`choices` below are the possible values from a GENCODE GTF file. If you are using a different GTF source, you may need to adjust this parameter.",
            choices: [
                "gene",
                "transcript",
                "exon",
                "CDS",
                "UTR",
                "start_codon",
                "stop_codon",
                "Selenocysteine",
            ],
        }
        protein: "Construct a protein database?"
        reference_fa_disk_size_gb: "Disk size (in GB) to allocate for downloading the reference FASTA file"
        gtf_disk_size_gb: "Disk size (in GB) to allocate for downloading the GTF file"
        kraken_fastas_disk_size_gb: "Disk size (in GB) to allocate for downloading the FASTA files"
    }

    input {
        String reference_fa_url
        String reference_fa_name
        String gtf_url
        String gtf_name
        Array[File] kraken_fastas = []
        Array[String] kraken_fasta_urls = []
        Array[String] kraken_libraries = [
            "archaea",
            "bacteria",
            "plasmid",
            "viral",
            "human",
            "fungi",
            "protozoa",
            "UniVec_Core",
        ]
        Array[String] coverage_feature_types = [
            "exon",
            "CDS",
            "UTR",
        ]
        Boolean protein = false
        Int reference_fa_disk_size_gb = 10
        Int gtf_disk_size_gb = 10
        Int kraken_fastas_disk_size_gb = 10
    }

    call util.download as reference_download { input:
        url = reference_fa_url,
        outfile_name = reference_fa_name,
        disk_size_gb = reference_fa_disk_size_gb,
    }
    call util.download as gtf_download { input:
        url = gtf_url,
        outfile_name = gtf_name,
        disk_size_gb = gtf_disk_size_gb,
    }

    scatter (feature_type in coverage_feature_types) {
        call util.make_coverage_regions_bed { input:
            gtf = gtf_download.downloaded_file,
            feature_type,
        }
    }

    scatter (url in kraken_fasta_urls) {
        call util.download as fastas_download { input:
            url,
            outfile_name = "tmp.fa.gz",
            disk_size_gb = kraken_fastas_disk_size_gb,
        }
    }

    scatter (lib in kraken_libraries) {
        call kraken2.download_library { input:
            library_name = lib,
            protein,
        }
    }

    if (
        (length(kraken_fastas) > 0)
        || (length(kraken_fasta_urls) > 0)
        || (length(kraken_libraries) > 0)
    ) {
        call kraken2.download_taxonomy { input: protein }
    }

    Array[File] custom_fastas = flatten([kraken_fastas, fastas_download.downloaded_file])
    if (length(custom_fastas) > 0) {
        call kraken2.create_library_from_fastas { input:
            fastas_gz = custom_fastas,
            protein,
        }
    }

    call kraken2.build_db as kraken_build_db { input:
        tarballs = flatten([
            [download_taxonomy.taxonomy],
            download_library.library,
            select_all([create_library_from_fastas.custom_library]),
        ]),
        protein,
    }

    output {
        File reference_fa = reference_download.downloaded_file
        File gtf = gtf_download.downloaded_file
        File kraken_db = kraken_build_db.built_db
        Array[File] coverage_beds = make_coverage_regions_bed.bed
    }
}
