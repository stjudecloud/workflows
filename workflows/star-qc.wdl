import "../tools/qc.wdl"
import "../tools/samtools.wdl"
import "../tools/picard.wdl"
import "../tools/fastqc.wdl"
import "../tools/qualimap.wdl"

workflow star_qc {
    File bam
    Int ncpu

    call samtools.quickcheck { input: bam=bam }
    call picard.mark_duplicates { input: bam=bam }
    call picard.validate_bam { input: bam=bam }
    call fastqc.fastqc { input: bam=bam, ncpu=ncpu }
    call qualimap.bamqc { input: bam=bam, ncpu=ncpu }
}
