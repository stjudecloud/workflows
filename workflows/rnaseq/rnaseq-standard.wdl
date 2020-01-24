## Description:
##
## This WDL workflow runs the STAR RNA-seq alignment workflow for St. Jude Cloud.
## The workflow takes an input BAM file and splits it into fastq files for each read in the pair. 
## The read pairs are then passed through STAR alignment to generate a BAM file. The BAM is run
## through several QC steps including FastQC and Qualimap. Quantification is done using htseq-count. 
## A final QC report is produced by MultiQC to generate a combined overview of the QC results
## for the sample.
##
## Inputs:
##
## reference_fasta - the genome for which to generate STAR reference files in FASTA format 
## gencode_gtf - the gene model file for the reference genome to use when generating STAR reference files 
## bam - input BAM file to realign
##
## LICENSING:
##
## MIT License
##
## Copyright 2019 St. Jude Children's Research Hospital
##
## Permission is hereby granted, free of charge, to any person obtaining a copy of this
## software and associated documentation files (the "Software"), to deal in the Software
## without restriction, including without limitation the rights to use, copy, modify, merge,
## publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons
## to whom the Software is furnished to do so, subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in all copies or
## substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
## BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
## NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
## DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
## OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

version 1.0

import "https://raw.githubusercontent.com/stjudecloud/workflows/master/workflows/bam-to-fastqs.wdl" as b2fq
import "https://raw.githubusercontent.com/stjudecloud/workflows/master/tools/star.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/master/tools/picard.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/master/tools/fastqc.wdl" as fqc
import "https://raw.githubusercontent.com/stjudecloud/workflows/master/tools/ngsderive.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/master/tools/qualimap.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/master/tools/htseq.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/master/tools/samtools.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/master/tools/md5sum.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/master/tools/multiqc.wdl" as mqc
import "https://raw.githubusercontent.com/stjudecloud/workflows/master/tools/qc.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/master/tools/util.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/master/tools/deeptools.wdl"

workflow rnaseq_standard {
    input {
        File gencode_gtf
        File input_bam
        File stardb_tar_gz
        String strand = ""
        String output_prefix = "out"
    }

    call util.get_read_groups { input: bam=input_bam }
    call util.prepare_read_groups_for_star { input: read_groups=get_read_groups.out }
    call b2fq.bam_to_fastqs { input: bam=input_bam }
    call star.alignment {
        input:
            read_one_fastqs=bam_to_fastqs.read1s,
            read_two_fastqs=bam_to_fastqs.read2s,
            stardb_tar_gz=stardb_tar_gz,
            output_prefix=output_prefix,
            read_groups=prepare_read_groups_for_star.out
    }
    call picard.sort as picard_sort { input: bam=alignment.star_bam }
    call samtools.index as samtools_index { input: bam=picard_sort.sorted_bam }
    call picard.validate_bam { input: bam=picard_sort.sorted_bam }
    call qc.parse_validate_bam { input: in=validate_bam.out }
    call fqc.fastqc { input: bam=picard_sort.sorted_bam }
    call ngsderive.infer_strand as ngsderive_strandedness { input: bam=picard_sort.sorted_bam, bai=samtools_index.bai, gtf=gencode_gtf }
    call qualimap.bamqc as qualimap_bamqc { input: bam=picard_sort.sorted_bam }
    call qualimap.rnaseq as qualimap_rnaseq { input: bam=picard_sort.sorted_bam, gencode_gtf=gencode_gtf }
    call htseq.count as htseq_count { input: bam=picard_sort.sorted_bam, gtf=gencode_gtf, strand=strand, inferred=ngsderive_strandedness.strandedness}
    call samtools.flagstat as samtools_flagstat { input: bam=picard_sort.sorted_bam }
    call md5sum.compute_checksum { input: infile=picard_sort.sorted_bam }
    call deeptools.bamCoverage as deeptools_bamCoverage { input: bam=picard_sort.sorted_bam, bai=samtools_index.bai }
    call mqc.multiqc {
        input:
            sorted_bam=picard_sort.sorted_bam,
            validate_sam_string=validate_bam.out,
            qualimap_bamqc=qualimap_bamqc.out_files,
            qualimap_rnaseq=qualimap_rnaseq.out_files,
            fastqc_files=fastqc.out_files,
            flagstat_file=samtools_flagstat.outfile,
            bigwig_file=deeptools_bamCoverage.bigwig,
            star_log=alignment.star_log
    }

    output {
        File bam = picard_sort.sorted_bam
        File bam_checksum = compute_checksum.outfile
        File bam_index = samtools_index.bai
        File bigwig = deeptools_bamCoverage.bigwig
        File gene_counts = htseq_count.out
        File flagstat = samtools_flagstat.outfile
        Array[File] fastqc_results = fastqc.out_files
        Array[File] qualimap_bamqc_results = qualimap_bamqc.out_files
        Array[File] qualimap_rnaseq_results = qualimap_rnaseq.out_files
        File multiqc_zip = multiqc.out
        File inferred_strandedness = ngsderive_strandedness.strandedness_file
    }
}
