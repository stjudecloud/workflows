version 1.1

import "../../tools/cellranger.wdl"
import "../../tools/md5sum.wdl"
import "../../tools/ngsderive.wdl"
import "../../tools/picard.wdl"
import "../../tools/samtools.wdl"
import "./10x-bam-to-fastqs.wdl" as bam_to_fastqs

workflow scrnaseq_standard {
    meta {
        description: "Align 10x Genomics FASTQ files to a reference genome and perform quantification."
        allowNestedInputs: true
        outputs: {
            harmonized_bam: "Aligned BAM file",
            bam_checksum: "Checksum of aligned BAM file",
            bam_index: "Index of aligned BAM file",
            qc: "Quality control metrics",
            barcodes: "Barcode information",
            features: "Feature information",
            matrix: "Gene expression matrix",
            filtered_gene_h5: "Filtered gene expression matrix",
            raw_gene_h5: "Raw gene expression matrix",
            raw_barcodes: "Raw barcode information",
            raw_features: "Raw feature information",
            raw_matrix: "Raw gene expression matrix",
            mol_info_h5: "Molecule information",
            web_summary: "HTML summary",
            inferred_strandedness: "Inferred strandedness",
        }
    }

    parameter_meta {
        bam: "Input BAM format file to quality check"
        gtf: "Gzipped GTF feature file"
        transcriptome_tar_gz: "Database of reference files for Cell Ranger. Can be downloaded from 10x Genomics."
        prefix: "Prefix for output files"
        validate_input: "Ensure input BAM is well-formed before beginning harmonization?"
        use_all_cores: "Use all cores for multi-core steps?"
        subsample_n_reads: "Only process a random sampling of `n` reads. <=`0` for processing entire input BAM."
    }

    input {
        File bam
        File gtf
        File transcriptome_tar_gz
        String prefix = basename(bam, ".bam")
        Boolean validate_input = true
        Boolean use_all_cores = false
        Int subsample_n_reads = -1
    }

    if (validate_input) {
        call picard.validate_bam as validate_input_bam { input:
            bam,
        }
    }

    if (subsample_n_reads > 0) {
        call samtools.subsample { input:
            bam,
            desired_reads = subsample_n_reads,
            use_all_cores,
        }
    }
    File selected_bam = select_first([subsample.sampled_bam, bam])

    call bam_to_fastqs.cell_ranger_bam_to_fastqs { input:
        bam = selected_bam,
        use_all_cores,
    }

    call cellranger.count { input:
        fastqs_tar_gz = cell_ranger_bam_to_fastqs.fastqs_archive,
        transcriptome_tar_gz,
        id = prefix,
        use_all_cores,
    }
    call picard.validate_bam { input: bam = count.bam }
    call ngsderive.strandedness { input:
        bam = count.bam,
        bam_index = count.bam_index,
        gene_model = gtf,
    }

    call md5sum.compute_checksum { input: file = count.bam }

    output {
        File harmonized_bam = count.bam
        File bam_checksum = compute_checksum.md5sum
        File bam_index = count.bam_index
        File qc = count.qc
        File barcodes = count.barcodes
        File features = count.features
        File matrix = count.matrix
        File filtered_gene_h5 = count.filtered_gene_h5
        File raw_gene_h5 = count.raw_gene_h5
        File raw_barcodes = count.raw_barcodes
        File raw_features = count.raw_features
        File raw_matrix = count.raw_matrix
        File mol_info_h5 = count.mol_info_h5
        File web_summary = count.web_summary
        File inferred_strandedness = strandedness.strandedness_file
    }
}
