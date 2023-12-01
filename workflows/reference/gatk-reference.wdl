# SPDX-License-Identifier: MIT
# Copyright St. Jude Children's Research Hospital
version 1.1

import "../../tools/picard.wdl"
import "../../tools/samtools.wdl"
import "../../tools/util.wdl"

workflow gatk_reference {
    meta {
        description: "Fetches reference files for GATK."
        outputs: {
            
        }
        allowNestedInputs: true
    }

    parameter_meta {

    }

    input {
        String reference_fa
        String reference_fa_name
        String reference_fa_md5
        String dbSNP_vcf_url
        String dbSNP_vcf_name
        String? dbSNP_vcf_index_url
        String? dbSNP_vcf_index_name
        Array[String] knownVCF_urls
        Array[String] knownVCF_names
        String? interval_list_url
        String? interval_list_name
    }

    call util.download as fasta_download {
        input:
            url = reference_fa,
            outfile_name = reference_fa_name,
            md5sum = reference_fa_md5
    }

    call samtools.faidx {
        input:
            fasta = fasta_download.downloaded_file
    }

    call picard.create_sequence_dictionary {
        input:
            fasta = faidx.indexed_fasta
    }

    call util.download as dbsnp {
        input:
            url = dbSNP_vcf_url,
            outfile_name = dbSNP_vcf_name
    }

    if (defined(dbSNP_vcf_index_url)) {
        call util.download as dbsnp_index {
            input:
                url = select_first([dbSNP_vcf_index_url, '']),
                outfile_name = select_first([dbSNP_vcf_index_name, ''])
        }
    }

    if (defined(interval_list_url)) {
        call util.download as intervals {
            input:
                url = select_first([interval_list_url, '']),
                outfile_name = select_first([interval_list_name, ''])
        }
    }

    scatter (pair in zip(knownVCF_urls, knownVCF_names)) {
        call util.download as knownVCF {
            input:
                url = pair.left,
                outfile_name = pair.right
        }
    }

    output {
        File fasta = faidx.indexed_fasta
        File fasta_index = faidx.fasta_index
        File fasta_dict = create_sequence_dictionary.dictionary
        File? dbSNP_vcf = dbsnp.downloaded_file
        File? dbSNP_vcf_index = dbsnp_index.downloaded_file
        File? interval_list = intervals.downloaded_file
        Array[File] knownVCFs = knownVCF.downloaded_file
    }
}