version 1.1

import "../../tools/kraken2.wdl"
import "../../tools/util.wdl"

workflow make_qc_reference {
    meta {
        description: "Downloads and creates all reference files needed to run the `quality_check` workflow"
        outputs: {
            reference_fa: "FASTA format reference file",
            gtf: "GTF feature file",
            exon_bed: "3 column BED file defining the regions of the exome. Derived from `gtf`.",
            CDS_bed: "3 column BED file defining the regions of the coding domain. Derived from `gtf`.",
            kraken_db: "A complete Kraken2 database"
        }
        allowNestedInputs: true
    }

    parameter_meta {
        reference_fa_url: "URL to retrieve the reference FASTA file from"
        reference_fa_name: "Name of output reference FASTA file"
        gtf_url: "URL to retrieve the reference GTF file from"
        gtf_name: "Name of output GTF file"
        kraken_fastas: "Array of gzipped FASTA files. Each sequence's ID must contain either an NCBI accession number or an explicit assignment of the taxonomy ID using `kraken:taxid`"
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
                "UniVec_Core"
            ]
        }
        kraken_fasta_urls: "URLs for any additional FASTA files in NCBI format to download and include in the Kraken2 database. This allows the addition of individual genomes (or other sequences) of interest."
        protein: "Construct a protein database?"
        kraken_fastas_disk_size_gb: "Disk size (in GB) to allocate for downloading the FASTA files"
        reference_fa_disk_size_gb: "Disk size (in GB) to allocate for downloading the reference FASTA file"
        gtf_disk_size_gb: "Disk size (in GB) to allocate for downloading the GTF file"
    }

    input {
        String reference_fa_url
        String reference_fa_name
        String gtf_url
        String gtf_name
        Array[File] kraken_fastas = []
        Array[String] kraken_libraries = [
            "archaea",
            "bacteria",
            "plasmid",
            "viral",
            "human",
            "fungi",
            "protozoa",
            "UniVec_Core"
        ]
        Array[String] kraken_fasta_urls = []
        Boolean protein = false
        Int kraken_fastas_disk_size_gb = 10
        Int reference_fa_disk_size_gb = 10
        Int gtf_disk_size_gb = 10
    }

    call util.download as reference_download { input:
        url=reference_fa_url,
        outfile_name=reference_fa_name,
        disk_size_gb=reference_fa_disk_size_gb,
    }
    call util.download as gtf_download { input:
        url=gtf_url,
        outfile_name=gtf_name,
        disk_size_gb=gtf_disk_size_gb,
    }

    call util.make_coverage_regions_beds { input:
        gtf=gtf_download.downloaded_file,
    }

    scatter (url in kraken_fasta_urls) {
        call util.download as fastas_download { input:
            url=url,
            outfile_name="tmp.fa.gz",
            disk_size_gb=kraken_fastas_disk_size_gb,
        }
    }

    call kraken2.download_taxonomy { input: protein=protein }

    scatter (lib in kraken_libraries) {
        call kraken2.download_library { input:
            library_name=lib,
            protein=protein,
        }
    }

    Array[File] custom_fastas = flatten([kraken_fastas, fastas_download.downloaded_file])
    Array[File] empty_array = []  # this structure is required by the WDL v1.1 spec
    if (custom_fastas != empty_array) {
        call kraken2.create_library_from_fastas { input:
            fastas_gz=custom_fastas,
            protein=protein,
        }
    }

    call kraken2.build_db as kraken_build_db { input:
        tarballs=flatten([
            [download_taxonomy.taxonomy],
            download_library.library,
            select_all([create_library_from_fastas.custom_library])
        ]),
        protein=protein,
    }

    output {
        File reference_fa = reference_download.downloaded_file
        File gtf = gtf_download.downloaded_file
        File exon_bed = make_coverage_regions_beds.exon_bed
        File CDS_bed = make_coverage_regions_beds.CDS_bed
        File kraken_db = kraken_build_db.built_db
    }
}
