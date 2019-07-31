import "./bam-to-fastqs.wdl" as b2fq
import "./star-db-build.wdl" as stardb_build
import "./star-alignment.wdl" as align
import "../tools/picard.wdl" 
import "../tools/fastqc.wdl"
import "../tools/rseqc.wdl"
import "../tools/qualimap.wdl"
import "../tools/htseq.wdl"
import "../tools/samtools.wdl"
import "../tools/md5sum.wdl"
import "../tools/multiqc.wdl"
import "../tools/qc.wdl"
import "../tools/util.wdl"

workflow start_to_finish {
    File reference_fasta 
    File gencode_gtf
    File refgene_bed
    File bam
    Int ncpu = 1

    call stardb_build.build_db {
        input:
            reference_fasta=reference_fasta,
            gencode_gtf=gencode_gtf,
            stardb_dir_name="stardb"
    }
    call util.get_read_groups { input: bam=bam }
    call util.prepare_read_groups_for_star { input: read_groups=get_read_groups.out }
    call b2fq.bam_to_fastqs { input: bam=bam }
    call align.star_alignment {
        input:
            read_one_fastqs=bam_to_fastqs.read1s,
            read_two_fastqs=bam_to_fastqs.read2s,
            stardb_dir=build_db.out_dir,
            output_prefix="out",
            read_groups=prepare_read_groups_for_star.out
    }
    call picard.mark_duplicates { input: bam=star_alignment.star_bam }
    call picard.validate_bam { input: bam=mark_duplicates.out }
    call qc.parse_validate_bam { input: in=validate_bam.out }
    call fastqc.fastqc { input: bam=mark_duplicates.out, ncpu=ncpu }   
    call rseqc.infer_experiment { input: bam=mark_duplicates.out, refgene_bed=refgene_bed}
    call qc.parse_infer_experiment { input: in=infer_experiment.out } 
    call qualimap.bamqc { input: bam=mark_duplicates.out, ncpu=ncpu }
    call qualimap.rnaseq { input: bam=mark_duplicates.out, gencode_gtf=gencode_gtf }
    call htseq.count { input: bam=mark_duplicates.out, gtf=gencode_gtf }
    call samtools.flagstat { input: bam=mark_duplicates.out }
    call samtools.index { input: bam=mark_duplicates.out }
    call md5sum.compute_checksum { input: infile=mark_duplicates.out }
    call multiqc.multiqc {
        input:
            star=star_alignment.star_bam,
            dups=mark_duplicates.out,
            validate_sam_string=validate_bam.out,
            qualimap_bamqc=bamqc.out_files,
            qualimap_rnaseq=rnaseq.out_files,
            fastqc_files=fastqc.out_files,
            flagstat_file=flagstat.flagstat
    }

}
