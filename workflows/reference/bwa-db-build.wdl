version 1.1

import "../../tools/bwa.wdl"
import "../../tools/util.wdl"

workflow bwa_db_build {
    meta {
        description: "Generates a set of genome reference files usable by the BWA aligner from an input reference file in FASTA format."
        outputs: {
            reference_fa: "FASTA format reference file used to generate `bwa_db_tar_gz`",
            bwa_db_tar_gz: "Gzipped tar archive of the BWA reference files. Files are at the root of the archive."
        }
        allowNestedInputs: true
    }

    parameter_meta {
        reference_fa_url: "URL to retrieve the reference FASTA file from."
        reference_fa_name: "Name of output reference FASTA file"
        reference_fa_md5: "Expected md5sum of reference FASTA file"
        reference_fa_disk_size_gb: "Disk size in GB to allocate for the reference FASTA file."
    }

    input {
        String reference_fa_url
        String reference_fa_name
        String? reference_fa_md5
        Int reference_fa_disk_size_gb = 10
    }

    call util.download as reference_download { input:
        url = reference_fa_url,
        outfile_name = reference_fa_name,
        disk_size_gb = reference_fa_disk_size_gb,
        md5sum = reference_fa_md5,
    }
    call bwa.build_bwa_db { input:
        reference_fasta = reference_download.downloaded_file,
    }

    output {
        File reference_fa = reference_download.downloaded_file
        File bwa_db_tar_gz = build_bwa_db.bwa_db_tar_gz
    }
}
