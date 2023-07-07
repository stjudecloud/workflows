## # Make QC Reference
##
## Download/create all the reference files needed for `quality-check-standard.wdl`.
## This includes: a reference FASTA, a GTF, the database used by Kraken2, and exonic/coding regions BEDs
## for use in restricting mosdepth coverage analysis to selected regions.
##
## ## LICENSING:
## 
## #### MIT License
##
## Copyright 2020-Present St. Jude Children's Research Hospital
##
## Permission is hereby granted, free of charge, to any person obtaining a copy of this
## software and associated documentation files (the "Software"), to deal in the Software
## without restriction, including without limitation the rights to use, copy, modify, merge,
## publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons
## to whom the Software is furnished to do so, subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in all copies or
## substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
## BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
## NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
## DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
## OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

version 1.1

import "../../tools/util.wdl"
import "../../tools/kraken.wdl"

workflow make_qc_reference {
    parameter_meta {
        kraken_fastas: "Array of gzipped FASTA files. Each sequence's ID must contain either an NCBI accession number or an explicit assignment of the taxonomy ID using `kraken:taxid`"
        kraken_libraries: {
            description: "List of kraken libraries to download"
            choices: [
                'archaea',
                'bacteria',
                'plasmid',
                'viral',
                'human',
                'fungi',
                'plant',
                'protozoa',
                'nt',
                'UniVec',
                'UniVec_Core'
            ]
        }
        kraken_fasta_urls: "URLs for any additional FASTA files in NCBI format to download and include in the Kraken2 database. This allows the addition of individual genomes (or other sequences) of interest."
        reference_fa_url: "URL to retrieve the reference FASTA file from"
        reference_fa_name: "Name of output reference FASTA file"
        gtf_url: "URL to retrieve the reference GTF file from"
        gtf_name: "Name of output GTF file"
        protein: "Construct a protein database?"
        max_retries: "Number of times to retry failed steps. Overrides task level defaults."
    }

    input {
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
        String reference_fa_url
        String reference_fa_name
        String gtf_url
        String gtf_name
        Boolean protein = false
        Int kraken_fastas_disk_size_gb = 10
        Int reference_fa_disk_size_gb = 10
        Int gtf_disk_size_gb = 10
        Int? max_retries
    }

    call util.download as reference_download { input:
        url=reference_fa_url,
        outfile_name=reference_fa_name,
        disk_size_gb=reference_fa_disk_size_gb,
        max_retries=max_retries
    }
    call util.download as gtf_download { input:
        url=gtf_url,
        outfile_name=gtf_name,
        disk_size_gb=gtf_disk_size_gb,
        max_retries=max_retries
    }

    call util.make_coverage_regions_beds { input:
        gtf=gtf_download.downloaded_file,
        max_retries=max_retries
    }

    scatter (url in kraken_fasta_urls) {
        call util.download as fastas_download { input:
            url=url,
            outfile_name="tmp.fa.gz",
            disk_size_gb=kraken_fastas_disk_size_gb,
            max_retries=max_retries
        }
    }

    call kraken.download_taxonomy { input: protein=protein, max_retries=max_retries }

    scatter (lib in kraken_libraries) {
        call kraken.download_library { input:
            library_name=lib,
            protein=protein,
            max_retries=max_retries
        }
    }

    Array[File] custom_fastas = flatten([kraken_fastas, fastas_download.downloaded_file])
    Array[File] empty_array = []  # this structure is required by the WDL v1 spec
    if (custom_fastas != empty_array) {
        call kraken.create_library_from_fastas { input:
            fastas=custom_fastas,
            protein=protein,
            max_retries=max_retries
        }
    }

    call kraken.build_db as kraken_build_db { input:
        tarballs=flatten([
            [download_taxonomy.taxonomy],
            download_library.library,
            select_all([create_library_from_fastas.custom_library])
        ]),
            protein=protein,
        max_retries=max_retries
    }

    output {
        File reference_fa = reference_download.downloaded_file
        File gtf = gtf_download.downloaded_file
        File exon_bed = make_coverage_regions_beds.exon_bed
        File CDS_bed = make_coverage_regions_beds.CDS_bed
        File kraken_db = kraken_build_db.built_db
    }
}
