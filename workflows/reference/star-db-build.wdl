# SPDX-License-Identifier: MIT
# Copyright St. Jude Children's Research Hospital
version 1.1

import "../../tools/star.wdl"
import "../../tools/util.wdl"

workflow star_db_build {
    meta {
        description: "Builds a database suitable for running the STAR alignment program"
        outputs: {
            reference_fa: "FASTA format reference file"
            gtf: "GTF feature file"
            star_db_tar_gz: "A gzipped TAR file containing the STAR reference files"
        }
        allowNestedInputs: true
    }

    parameter_meta {
        reference_fa_url: "URL to retrieve the reference FASTA file from"
        reference_fa_name: "Name of output reference FASTA file"
        gtf_url: "URL to retrieve the reference GTF file from"
        gtf_name: "Name of output GTF file"
        reference_fa_md5: "Expected md5sum of reference FASTA file"
        gtf_md5: "Expected md5sum of GTF file"
        reference_fa_disk_size_gb: "Disk space to allocate the FASTA download task"
        gtf_disk_size_gb: "Disk space to allocate the GTF download task"
    }

    input {
        String reference_fa_url
        String reference_fa_name
        String gtf_url
        String gtf_name
        String? reference_fa_md5
        String? gtf_md5
        Int reference_fa_disk_size_gb = 10
        Int gtf_disk_size_gb = 10
    }

    call util.download as reference_download { input:
        url=reference_fa_url,
        outfile_name=reference_fa_name,
        md5sum=reference_fa_md5,
        disk_size_gb=reference_fa_disk_size_gb,
    }
    call util.download as gtf_download { input:
        url=gtf_url,
        outfile_name=gtf_name,
        md5sum=gtf_md5,
        disk_size_gb=gtf_disk_size_gb,
    }
    call star.build_star_db { input:
        reference_fasta=reference_download.downloaded_file,
        gtf=gtf_download.downloaded_file,
    }

    output {
      File reference_fa = reference_download.downloaded_file
      File gtf = gtf_download.downloaded_file
      File star_db_tar_gz = build_star_db.star_db
    }
}
