version 1.2

import "../../data_structures/read_group.wdl" as read_group_ds
import "../../tools/bwa.wdl" as bwa
import "../../tools/bwamem2.wdl" as bwamem2
import "../../tools/minimap2.wdl" as minimap2
import "../../tools/samtools.wdl" as samtools
import "../../tools/vg.wdl" as vg
import "../general/alignment-post.wdl" as alignment_post

workflow align {
    meta {
        name: "Aligner comparison workflow"
        description: "Align sequencing reads to a reference genome with a variety of aligners"
        category: "Utility"
        outputs: {
            bwamem2_bam: "BAM file output from BWA-MEM2 alignment",
            bwamem2_bam_index: "BAM index file for BWA-MEM2 alignments",
            vg_giraffe_bam: "BAM file output from vg giraffe alignment",
            vg_giraffe_bam_index: "BAM index file for vg giraffe alignments",
            minimap2_bam: "BAM file output from minimap2 alignment",
            minimap2_bam_index: "BAM index file for minimap2 alignments",
            bwa_bam: "BAM file output from BWA alignment",
            bwa_bam_index: "BAM index file for BWA alignments",
        }
        allowNestedInputs: true
    }

    parameter_meta {
        read_one_fastq_gz: "Input gzipped FASTQ read one file to align"
        bwamem2_reference: "The BWA-MEM2 index file for the reference genome"
        giraffe_gbz_graph: "The vg GBZ graph file for the reference genome"
        giraffe_minimizer_index: "The vg minimizer index file for the reference genome"
        giraffe_zipcode_name: "The vg zipcode name file for the reference genome"
        giraffe_distance_index: "The vg distance index file for the reference genome"
        minimap2_reference: "The minimap2 index file for the reference genome"
        bwamem_reference: "The BWA index tar.gz file for the reference genome"
        read_group: "The read group to be included in the SAM header"
        read_two_fastq_gz: "Input gzipped FASTQ read two file to align (for paired-end data)"
        output_prefix: "Output file prefix for aligned reads"
        threads: "Number of threads to use for alignment"
        modify_disk_size_gb: "Additional disk space to allocate (in GB)"
    }

    input {
        File read_one_fastq_gz
        File bwamem2_reference
        File giraffe_gbz_graph
        File giraffe_minimizer_index
        File giraffe_zipcode_name
        File giraffe_distance_index
        File minimap2_reference
        File bwamem_reference
        ReadGroup read_group
        File? read_two_fastq_gz
        String output_prefix = "aligned_output"
        Int threads = 30
        Int modify_disk_size_gb = 0
    }

    call read_group_ds.read_group_to_string as read_group_string {
        read_group,
        format_as_sam_record = true,
    }

    call bwamem2.align as bwamem2 {
        read_one_fastq_gz,
        read_two_fastq_gz,
        reference_index = bwamem2_reference,
        prefix = "~{output_prefix}_bwamem2",
        threads,
        modify_disk_size_gb,
        read_group = read_group_string.validated_read_group,
    }

    call alignment_post.alignment_post as bwamem2_post {
        bam = bwamem2.alignments,
        mark_duplicates = true,
    }

    call vg.giraffe as vg_giraffe {
        read_one_fastq_gz,
        read_two_fastq_gz,
        gbz_graph = giraffe_gbz_graph,
        minimizer_index = giraffe_minimizer_index,
        zipcode_name = giraffe_zipcode_name,
        distance_index = giraffe_distance_index,
        output_name = "~{output_prefix}_vg.bam",
        output_format = "BAM",
        threads,
        modify_disk_size_gb,
    }

    call read_group_ds.read_group_to_array {
        read_group,
    }

    call samtools.addreplacerg {
        bam = vg_giraffe.alignments,
        orphan_only = false,
        read_group_line = read_group_to_array.converted_read_group,
    }

    call alignment_post.alignment_post as giraffe_post {
        bam = addreplacerg.tagged_bam,
        mark_duplicates = true,
    }

    call minimap2.align as minimap2 {
        read_one_fastq_gz,
        read_two_fastq_gz,
        reference_index = minimap2_reference,
        output_name = "~{output_prefix}_minimap2.bam",
        threads,
        modify_disk_size_gb,
        read_group = read_group_string.validated_read_group,
    }

    call alignment_post.alignment_post as minimap2_post {
        bam = minimap2.alignments,
        mark_duplicates = true,
    }

    call bwa.bwa_mem as bwa_mem {
        read_one_fastq_gz,
        read_two_fastq_gz,
        bwa_db_tar_gz = bwamem_reference,
        prefix = "~{output_prefix}_bwa",
        ncpu = threads,
        modify_disk_size_gb,
        read_group = read_group_string.validated_read_group,
    }

    call alignment_post.alignment_post as bwa_mem_post {
        bam = bwa_mem.bam,
        mark_duplicates = true,
    }

    output {
        File bwamem2_bam = bwamem2_post.processed_bam
        File bwamem2_bam_index = bwamem2_post.bam_index
        File vg_giraffe_bam = giraffe_post.processed_bam
        File vg_giraffe_bam_index = giraffe_post.bam_index
        File minimap2_bam = minimap2_post.processed_bam
        File minimap2_bam_index = minimap2_post.bam_index
        File bwa_bam = bwa_mem_post.processed_bam
        File bwa_bam_index = bwa_mem_post.bam_index
    }
}
