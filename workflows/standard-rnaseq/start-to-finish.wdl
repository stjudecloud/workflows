## Description: 
##
## This WDL workflow runs the STAR RNA-seq alignment workflow for St. Jude Cloud. The workflow takes an input BAM file and splits it into fastq files for each read in the pair. 
## The read pairs are then passed through STAR alignment to generate a BAM file. The BAM is run through several QC steps including FastQC and Qualimap.
## Quantification is done using htseq-count. A final QC report is produced by MultiQC to generate a combined overview of the QC results for the sample.  
##
## Inputs: 
##
## reference_fasta - the genome for which to generate STAR reference files in FASTA format 
## gencode_gtf - the gene model file for the reference genome to use when generating STAR reference files 
## refgene_bed - refgene variants in BED format
## bam - input BAM file to realign
## ncpu (optional) - the number of CPUs to use when running steps that support multithreading 
##
## LICENSING :
## MIT License
##
## Copyright 2019 St. Jude Children's Research Hospital
##
## Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
##
##The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
##
##THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import "https://raw.githubusercontent.com/stjudecloud/workflows/jobin/rnaseq_v2_azure/workflows/bam-to-fastqs.wdl" as b2fq
import "https://raw.githubusercontent.com/stjudecloud/workflows/jobin/rnaseq_v2_azure/workflows/star-db-build.wdl" as stardb_build
import "https://raw.githubusercontent.com/stjudecloud/workflows/jobin/rnaseq_v2_azure/workflows/star-alignment.wdl" as align
import "https://raw.githubusercontent.com/stjudecloud/workflows/jobin/rnaseq_v2_azure/tools/picard.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/jobin/rnaseq_v2_azure/tools/fastqc.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/jobin/rnaseq_v2_azure/tools/rseqc.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/jobin/rnaseq_v2_azure/tools/qualimap.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/jobin/rnaseq_v2_azure/tools/htseq.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/jobin/rnaseq_v2_azure/tools/samtools.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/jobin/rnaseq_v2_azure/tools/md5sum.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/jobin/rnaseq_v2_azure/tools/multiqc.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/jobin/rnaseq_v2_azure/tools/qc.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/jobin/rnaseq_v2_azure/tools/util.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/jobin/rnaseq_v2_azure/tools/deeptools.wdl"

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
            stardb_zip=build_db.stardb_zip,
            output_prefix="out",
            read_groups=prepare_read_groups_for_star.out
    }
    call picard.validate_bam { input: bam=star_alignment.star_bam }
    call qc.parse_validate_bam { input: in=validate_bam.out }
    call fastqc.fastqc { input: bam=star_alignment.star_bam, ncpu=ncpu }   
    call rseqc.infer_experiment { input: bam=star_alignment.star_bam, refgene_bed=refgene_bed}
    call qc.parse_infer_experiment { input: in=infer_experiment.out } 
    call qualimap.bamqc { input: bam=star_alignment.star_bam, ncpu=ncpu }
    call qualimap.rnaseq { input: bam=star_alignment.star_bam, gencode_gtf=gencode_gtf }
    call htseq.count { input: bam=star_alignment.star_bam, gtf=gencode_gtf }
    call samtools.flagstat { input: bam=star_alignment.star_bam }
    call samtools.index { input: bam=star_alignment.star_bam }
    call md5sum.compute_checksum { input: infile=star_alignment.star_bam } 
    call deeptools.bamCoverage { input: bam=star_alignment.star_bam, bai=index.bai }
    call multiqc.multiqc {
        input:
            star=star_alignment.star_bam,
            validate_sam_string=validate_bam.out,
            qualimap_bamqc=bamqc.out_files,
            qualimap_rnaseq=rnaseq.out_files,
            fastqc_files=fastqc.out_files,
            flagstat_file=flagstat.flagstat
    }
    output {
        star_alignment.star_bam
        index.bai
        count.out
        flagstat.flagstat
        multiqc.out
    }
}
