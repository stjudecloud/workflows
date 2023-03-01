## Cell Ranger
##
## This WDL tool wrap the [10x Genomics Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) tool.
## Cell Ranger is a tool for handling scRNA-Seq data.

version 1.0

task count {
    input {
        String id
        File transcriptome_tar_gz
        File fastqs_tar_gz
        Int ncpu = 8
        Int memory_gb = 16
        String jobmode = "local"
        Int max_retries = 1
        Boolean detect_nproc = false
    }

    Float fastq_size = size(fastqs_tar_gz, "GiB")
    Int disk_size = ceil((fastq_size * 2) + 10)
    String parsed_detect_nproc = if detect_nproc then "true" else ""

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if [ -n ~{parsed_detect_nproc} ]
        then
            n_cores=$(nproc)
        fi

        mkdir transcriptome_dir
        tar zxf ~{transcriptome_tar_gz} -C transcriptome_dir --strip-components 1
   
        mkdir fastqs
        tar zxf ~{fastqs_tar_gz} -C fastqs

        files=(fastqs/*.fastq.gz)
        # sample parameter to cellranger count must match the sample prefix contained in the fastq file.
        # So we infer it here by manipulating the file name.
        # expected sample name extension comes from:
        # https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq2-v2-20-software-guide-15051736-03.pdf
        sample_id="$(basename "${files[0]}" | sed -E 's/_S[1-9]_L[0-9]{3}_[I,R][1,2]_001.fastq.gz$//')"

        cellranger count \
            --id ~{id} \
            --transcriptome transcriptome_dir \
            --fastqs fastqs \
            --sample "${sample_id}" \
            --jobmode ~{jobmode} \
            --localcores "$n_cores" \
            --localmem ~{memory_gb} \
            --disable-ui

        ls > /dev/stderr
    >>>

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        docker: "ghcr.io/stjudecloud/cellranger:1.1.0"
        maxRetries: max_retries
        cpu: ncpu
    }

    output {
        File bam = glob("*/outs/possorted_genome_bam.bam")[0]
        File bam_index = glob("*/outs/possorted_genome_bam.bam.bai")[0]
        File qc = glob("*/outs/metrics_summary.csv")[0]
        File barcodes = glob("*/outs/filtered_feature_bc_matrix/barcodes.tsv.gz")[0]
        File features = glob("*/outs/filtered_feature_bc_matrix/features.tsv.gz")[0]
        File matrix = glob("*/outs/filtered_feature_bc_matrix/matrix.mtx.gz")[0]
        File filtered_gene_h5 = glob("*/outs/filtered_feature_bc_matrix.h5")[0]
        File raw_gene_h5 = glob("*/outs/raw_feature_bc_matrix.h5")[0]
        File raw_barcodes = glob("*/outs/raw_feature_bc_matrix/barcodes.tsv.gz")[0]
        File raw_features = glob("*/outs/raw_feature_bc_matrix/features.tsv.gz")[0]
        File raw_matrix = glob("*/outs/raw_feature_bc_matrix/matrix.mtx.gz")[0]
        File mol_info_h5 = glob("*/outs/molecule_info.h5")[0]
        File web_summary = glob("*/outs/web_summary.html")[0]
        File loupe = glob("*/outs/cloupe.cloupe" )[0]
    }

    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
        description: "This WDL tool runs Cell Ranger count to generate an aligned BAM and feature counts from scRNA-Seq data."
    }

    parameter_meta {
        id: "A unique run ID"
        fastqs_tar_gz: "Path to the fastq folder archive in .tar.gz format"
        transcriptome_tar_gz: "Path to Cell Ranger-compatible transcriptome reference in .tar.gz format"
        sample_id: "Sample name as used by cellranger mkfastq"
    }
}

task bamtofastq {
    input {
        File bam
        Int ncpu = 4
        Int memory_gb = 8
        Boolean cellranger11 = false
        Boolean longranger20 = false
        Boolean gemcode = false
        Boolean detect_nproc = false
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil((bam_size * 2) + 10)

    String data_arg = if (cellranger11) then "--cr11"
                        else if (longranger20) then "--lr10"
                        else if (gemcode) then "--gemcode"
                        else ""
    String parsed_detect_nproc = if detect_nproc then "true" else ""

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if [ -n ~{parsed_detect_nproc} ]
        then
            n_cores=$(nproc)
        fi
        
        cellranger bamtofastq --nthreads "$n_cores" ~{data_arg} ~{bam} fastqs
        cd fastqs/*/
        tar -zcf archive.tar.gz ./*.fastq.gz
    >>>

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        docker: "ghcr.io/stjudecloud/cellranger:1.1.0"
        maxRetries: max_retries
        cpu: ncpu
    }

    output {
        Array[File] fastqs = glob("fastqs/*/*fastq.gz")
        File fastqs_archive = glob("fastqs/*/*.tar.gz")[0]
        Array[File] read1 = glob("fastqs/*/*R1*.fastq.gz")
        Array[File] read2 = glob("fastqs/*/*R2*.fastq.gz")
    }

    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
        description: "This WDL tool runs the 10x bamtofastq tool to convert Cell Ranger generated BAM files back to fastq files"

    }

    parameter_meta {
        bam: "Input BAM to convert to Cell Ranger compatible fastqs"
        cellranger11: "Convert a BAM produced by Cell Ranger 1.0-1.1"
        longranger20: "Convert a BAM produced by Longranger 2.0"
        gemcode: "Convert a BAM produced from GemCode data (Longranger 1.0 - 1.3)"
    }
}
