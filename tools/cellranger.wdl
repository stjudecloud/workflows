## Cell Ranger
##
## This WDL file wrap the [10x Genomics Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) tool.
## Cell Ranger is a tool for handling scRNA-Seq data.
#@ except: LineWidth

version 1.1

task count {
    meta {
        description: "This WDL task runs Cell Ranger count to generate an aligned BAM and feature counts from scRNA-Seq data."
        outputs: {
            bam: "Aligned BAM file",
            bam_index: "BAM index file",
            qc: "Quality control metrics in CSV format",
            barcodes: "Barcodes in gzipped TSV format",
            features: "Filtered features in gzipped TSV format",
            matrix: "Filtered matrix of features",
            filtered_gene_h5: "Filtered gene matrix in H5 format",
            raw_gene_h5: "Raw gene matrix in H5 format",
            raw_barcodes: "Raw barcodes in gzipped TSV format",
            raw_features: "Raw features in gzipped TSV format",
            raw_matrix: "Raw matrix of features",
            mol_info_h5: "Molecule information in H5 format",
            web_summary: "HTML summary of the run",
            cloupe: "Cloupe file for visualization",
        }
    }

    parameter_meta {
        fastqs_tar_gz: "Path to the FASTQ folder archive in .tar.gz format"
        transcriptome_tar_gz: "Path to Cell Ranger-compatible transcriptome reference in .tar.gz format"
        id: "A unique run ID"
        use_all_cores: "Use all cores? Recommended for cloud environments."
        ncpu: "Number of cores to allocate for task"
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File fastqs_tar_gz
        File transcriptome_tar_gz
        String id
        Boolean use_all_cores = false
        Int ncpu = 1
        Int memory_gb = 16
        Int modify_disk_size_gb = 0
    }

    Float fastq_size = size(fastqs_tar_gz, "GiB")
    Float transcriptome_size = size(transcriptome_tar_gz, "GiB")
    Int disk_size_gb = (
        ceil((fastq_size + transcriptome_size) * 2) + 10 + modify_disk_size_gb
    )

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi

        mkdir transcriptome_dir
        tar -xzf ~{transcriptome_tar_gz} \
            -C transcriptome_dir \
            --strip-components 1 \
            --no-same-owner

        mkdir fastqs
        tar -xzf ~{fastqs_tar_gz} -C fastqs --no-same-owner

        files=(fastqs/*.fastq.gz)
        # sample parameter to cellranger count must match
        # the sample prefix contained in the FASTQ file.
        # So we infer it here by manipulating the file name.
        # expected sample name extension comes from:
        # https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq2-v2-20-software-guide-15051736-03.pdf
        sample_id="$(basename "${files[0]}" | sed -E 's/_S[1-9]_L[0-9]{3}_[I,R][1,2]_001.fastq.gz$//')"

        cellranger count \
            --id ~{id} \
            --transcriptome transcriptome_dir \
            --fastqs fastqs \
            --sample "${sample_id}" \
            --jobmode local \
            --localcores "$n_cores" \
            --localmem ~{memory_gb} \
            --disable-ui
    >>>

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
        File cloupe = glob("*/outs/cloupe.cloupe" )[0]
    }

    runtime {
        cpu: ncpu
        memory: "~{memory_gb} GB"
        disks: "~{disk_size_gb} GB"
        container: "ghcr.io/stjudecloud/cellranger:1.1.1"
        maxRetries: 1
    }
}

task bamtofastq {
    meta {
        description: "This WDL task runs the 10x bamtofastq tool to convert Cell Ranger generated BAM files back to FASTQ files"
        outputs: {
            fastqs: "FASTQ files",
            fastqs_archive: "FASTQ files in a tarball",
            read_one_fastq_gz: "Read 1 FASTQ files",
            read_two_fastq_gz: "Read 2 FASTQ files",
        }
    }

    parameter_meta {
        bam: "Input BAM to convert to Cell Ranger compatible fastqs"
        cellranger11: "Convert a BAM produced by Cell Ranger 1.0-1.1"
        longranger20: "Convert a BAM produced by Longranger 2.0"
        gemcode: "Convert a BAM produced from GemCode data (Longranger 1.0 - 1.3)"
        use_all_cores: "Use all cores? Recommended for cloud environments."
        ncpu: "Number of cores to allocate for task"
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        Boolean cellranger11 = false
        Boolean longranger20 = false
        Boolean gemcode = false
        Boolean use_all_cores = false
        Int ncpu = 1
        Int memory_gb = 40
        Int modify_disk_size_gb = 0
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size * 2) + 10 + modify_disk_size_gb

    String data_arg = if (cellranger11) then "--cr11"
                        else if (longranger20) then "--lr10"
                        else if (gemcode) then "--gemcode"
                        else ""

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi

        cellranger bamtofastq --nthreads "$n_cores" ~{data_arg} ~{bam} fastqs

        cd fastqs/*/
        tar -czf archive.tar.gz ./*.fastq.gz
    >>>

    output {
        Array[File] fastqs = glob("fastqs/*/*fastq.gz")
        File fastqs_archive = glob("fastqs/*/*.tar.gz")[0]
        Array[File] read_one_fastq_gz = glob("fastqs/*/*R1*.fastq.gz")
        Array[File] read_two_fastq_gz = glob("fastqs/*/*R2*.fastq.gz")
    }

    runtime {
        cpu: ncpu
        memory: "~{memory_gb} GB"
        disks: "~{disk_size_gb} GB"
        container: "ghcr.io/stjudecloud/cellranger:1.1.1"
        maxRetries: 1
    }
}
