import "../tools/qc.wdl"

workflow star_qc {
    File bam
    Int ncpu

    call qc.quickcheck { input: bam=bam }
    call qc.mark_duplicates { input: bam=bam }
    call qc.validate_bam { input: bam=bam }
    call qc.fastqc { input: bam=bam, ncpu=ncpu }
    call qc.qualimap { input: bam=bam, ncpu=ncpu }
}