version 1.1

import "../../tools/bowtie2.wdl"
import "../../tools/samtools.wdl"
import "../../tools/util.wdl"

workflow make_hic_reference {
    meta {
        description: "Downloads and creates all reference files needed to run the `hic` workflow"
        outputs: {
            reference_fasta: "FASTA format reference file used to generate `bowtie2_db_tar_gz`",
            reference_fasta_index: "FASTA index file for `reference_fasta`",
            reference_chromsizes: "Chromosome sizes file for `reference_fasta`",
            genome_fragment: "BED file(s) with restriction fragments",
            exclude_list: "BED file with regions to exclude from analysis",
            bowtie2_db_tar_gz: "Gzipped tar archive of the Bowtie2 reference files. Files are at the root of the archive.",
            bowtie2_index_files: "Bowtie2 index files",
        }
        allowNestedInputs: true
    }

    parameter_meta {
        reference_fa_url: "URL to retrieve the reference FASTA file from"
        reference_fa_name: "Name of output reference FASTA file"
        exclude_list_url: "URL from which to retrieve the exclude list file"
        exclude_list_name: "Name of output exclude list file"
        restriction_sites: {
            description: "List of restriction sites for which to extract restriction fragments",
            help: "This uses HiC-Pro's `digest_genome.py` script to generate restriction fragments. Each site should be specified as a string of nucleotides. A carat (`^`) marks the cut site.",
            external_help: "https://github.com/nservant/HiC-Pro/blob/master/doc/UTILS.md#digest_genomepy-or-how-can-i-generate-the-list-of-restriction-fragments-after-genome-digestion-",
        }
        restriction_sites_names: "Names for the restriction sites to use in output files"
    }

    input {
        Array[String] restriction_sites
        Array[String] restriction_sites_names
        String reference_fa_url
        String reference_fa_name
        String exclude_list_url
        String exclude_list_name
    }

    call parse_input { input:
        restriction_sites,
        restriction_sites_names,
    }

    call util.download as reference_download { input:
        url = reference_fa_url,
        outfile_name = reference_fa_name,
        disk_size_gb = 10,
    }

    call util.download as exclude_list_download { input:
        url = exclude_list_url,
        outfile_name = exclude_list_name,
        disk_size_gb = 1,
    }

    call bowtie2.build { input:
        reference = reference_download.downloaded_file,
        prefix = basename(reference_fa_name, ".fa.gz"),
    }

    call samtools.faidx { input:
        fasta = reference_download.downloaded_file
    }

    call chromsizes { input:
        fasta_index = faidx.fasta_index
    }

    scatter (site in zip(restriction_sites, restriction_sites_names)) {
        call fragment_file after parse_input { input:
            reference_fasta = reference_download.downloaded_file,
            restriction_site = site.left,
            output_name = basename(reference_fa_name, ".gz") + "." + site.right + ".bed"
        }
    }

    output {
        File reference_fasta = reference_download.downloaded_file
        File reference_fasta_index = faidx.fasta_index
        File reference_chromsizes = chromsizes.chromsizes
        Array[File]? genome_fragment = fragment_file.fragment_file
        File exclude_list = exclude_list_download.downloaded_file
        File bowtie2_db_tar_gz = build.bowtie_db_tar_gz
        Array[File] bowtie2_index_files = build.index_files
    }
}

task chromsizes {
    meta {
        description: "Create a chromsizes file from a fasta index file"
        outputs: {
            chromsizes: "TSV with chromosome names and sizes"
        }
    }

    parameter_meta {
        fasta_index: "Fasta index file from which to create chromsizes file"
        output_name: "Name of output file containing chromosome sizes"
    }

    input {
        File fasta_index
        String output_name = basename(fasta_index, ".fai") + ".tab"
    }

    command <<<
        cut -f 1,2 ~{fasta_index} > ~{output_name}
    >>>

    output {
        File chromsizes = output_name
    }

    runtime {
        cpu: 1
        memory: "4 GB"
        container: "ubuntu:22.04"
        maxRetries: 1
    }
}

task fragment_file {
    meta {
        description: "Create a fragment file from a reference FASTA file and a list of restriction site"
        outputs: {
            fragment_file: "BED file with restriction fragments"
        }
    }

    parameter_meta {
        reference_fasta: "Reference FASTA file"
        restriction_site: "Restriction site"
        output_name: "Name of output fragment file"
    }

    input {
        File reference_fasta
        String restriction_site
        String output_name
    }

    String base = basename(reference_fasta, ".gz")

    command <<<
        gunzip -c ~{reference_fasta} > ~{base} \
           || ln -sf ~{reference_fasta} ~{base}

        /HiC-Pro_3.0.0/bin/utils/digest_genome.py \
            -r ~{restriction_site} \
            -o ~{output_name} \
            ~{base}
    >>>

    output {
        File fragment_file = output_name
    }

    runtime {
        cpu: 1
        memory: "4 GB"
        container: "nservant/hicpro:3.0.0"
        maxRetries: 1
    }
}

task parse_input {
    meta {
        description: "Parse parameters for the `make_hic_reference` workflow"
        outputs: {
            check: "Dummy output to indicate success and to enable call-caching"
        }
    }

    parameter_meta {
        restriction_sites: {
            description: "List of restriction sites for which to extract restriction fragments",
            help: "This uses HiC-Pro's `digest_genome.py` script to generate restriction fragments. Each site should be specified as a string of nucleotides. A carat (`^`) marks the cut site.",
            external_help: "https://github.com/nservant/HiC-Pro/blob/master/doc/UTILS.md#digest_genomepy-or-how-can-i-generate-the-list-of-restriction-fragments-after-genome-digestion-",
        }
        restriction_sites_names: "Names for the restriction sites to use in output files"
    }

    input {
        Array[String] restriction_sites
        Array[String] restriction_sites_names
    }

    command <<<
        if [ ~{length(restriction_sites)} -ne ~{length(restriction_sites_names)} ]
        then
            >&2 echo -n "Length of restriction_sites "
            >&2 echo "and restriction_sites_names must be equal"
            exit 1
        fi
    >>>

    output {
        String check = "passed"
    }

    runtime {
        memory: "4 GB"
        disks: "10 GB"
        container: "ghcr.io/stjudecloud/util:1.3.0"
        maxRetries: 1
    }
}
