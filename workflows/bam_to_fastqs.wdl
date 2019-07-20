import "../tools/samtools.wdl"
import "../tools/picard.wdl"
import "../tools/fq.wdl"

workflow bam_to_fastqs {
    File bam

    call samtools.quickcheck { input: bam=bam }
    call samtools.split { input: bam=bam }
    scatter (split_bam in split.split_bams) {
        call picard.bam_to_fastq { input: bam=split_bam }
    }
    scatter (reads in zip(bam_to_fastq.read1, bam_to_fastq.read2)) {
        call fq.fqlint { input: read1=reads.left, read2=reads.right}
    }
}
