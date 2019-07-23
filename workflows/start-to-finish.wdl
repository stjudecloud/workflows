import "./bam-to-fastqs.wdl" as b2fq
import "./star-db-build.wdl" as stardb_build
import "./star-alignment.wdl" as align

workflow start_to_finish {
    File reference_fasta 
    File gencode_gtf
    File bam

    call stardb_build.build_db {
        input:
            reference_fasta=reference_fasta,
            gencode_gtf=gencode_gtf,
            stardb_dir_name="stardb"
    }
    call b2fq.bam_to_fastqs { input: bam=bam }
    call align.star_alignment {
        input:
            read_one_fastqs=bam_to_fastqs.read1s,
            read_two_fastqs=bam_to_fastqs.read2s,
            stardb_dir=build_db.out_dir,
            output_prefix="out"
    }
}