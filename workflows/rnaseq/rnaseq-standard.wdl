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

import "https://raw.githubusercontent.com/stjudecloud/workflows/master/workflows/general/bam-to-fastqs.wdl" as b2fq
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
        Int max_retries = 1
    }

    call parse_input { input: input_bam=input_bam, input_strand=strand }

    call util.get_read_groups { input: bam=parse_input.bam_after_input_validated, max_retries=max_retries }
    call util.prepare_read_groups_for_star { input: read_groups=get_read_groups.out, max_retries=max_retries }
    call b2fq.bam_to_fastqs { input: bam=parse_input.bam_after_input_validated, max_retries=max_retries }
    call star.alignment {
        input:
            read_one_fastqs=bam_to_fastqs.read1s,
            read_two_fastqs=bam_to_fastqs.read2s,
            stardb_tar_gz=stardb_tar_gz,
            output_prefix=output_prefix,
            read_groups=prepare_read_groups_for_star.out,
            max_retries=max_retries
    }
    call picard.sort as picard_sort { input: bam=alignment.star_bam, max_retries=max_retries }
    call samtools.index as samtools_index { input: bam=picard_sort.sorted_bam, max_retries=max_retries }
    call picard.validate_bam { input: bam=picard_sort.sorted_bam, max_retries=max_retries }
    call qc.parse_validate_bam { input: in=validate_bam.out, max_retries=max_retries }
    call fqc.fastqc { input: bam=picard_sort.sorted_bam, max_retries=max_retries }
    call ngsderive.infer_strand as ngsderive_strandedness { input: bam=picard_sort.sorted_bam, bai=samtools_index.bai, gtf=gencode_gtf, max_retries=max_retries }
    call qualimap.bamqc as qualimap_bamqc { input: bam=picard_sort.sorted_bam, max_retries=max_retries }
    call qualimap.rnaseq as qualimap_rnaseq { input: bam=picard_sort.sorted_bam, gencode_gtf=gencode_gtf, strand=strand, inferred=ngsderive_strandedness.strandedness, max_retries=max_retries }
    call htseq.count as htseq_count { input: bam=picard_sort.sorted_bam, gtf=gencode_gtf, strand=strand, inferred=ngsderive_strandedness.strandedness, max_retries=max_retries }
    call samtools.flagstat as samtools_flagstat { input: bam=picard_sort.sorted_bam, max_retries=max_retries }
    call md5sum.compute_checksum { input: infile=picard_sort.sorted_bam, max_retries=max_retries }
    call deeptools.bamCoverage as deeptools_bamCoverage { input: bam=picard_sort.sorted_bam, bai=samtools_index.bai, max_retries=max_retries }
    call mqc.multiqc {
        input:
            sorted_bam=picard_sort.sorted_bam,
            validate_sam_string=validate_bam.out,
            qualimap_bamqc=qualimap_bamqc.results,
            qualimap_rnaseq=qualimap_rnaseq.results,
            fastqc_files=fastqc.out_files,
            flagstat_file=samtools_flagstat.outfile,
            bigwig_file=deeptools_bamCoverage.bigwig,
            star_log=alignment.star_log, 
            max_retries=max_retries
    }

    output {
        File bam = picard_sort.sorted_bam
        File bam_checksum = compute_checksum.outfile
        File bam_index = samtools_index.bai
        File bigwig = deeptools_bamCoverage.bigwig
        File gene_counts = htseq_count.out
        File flagstat = samtools_flagstat.outfile
        Array[File] fastqc_results = fastqc.out_files
        File qualimap_bamqc_results = qualimap_bamqc.results
        File qualimap_rnaseq_results = qualimap_rnaseq.results
        File multiqc_zip = multiqc.out
        File inferred_strandedness = ngsderive_strandedness.strandedness_file
    }
}

task parse_input {
    input {
        File input_bam
        String input_strand
    }

    command {
        if [ -n "~{input_strand}" ] && [ "~{input_strand}" != "reverse" ] && [ "~{input_strand}" != "yes" ] && [ "~{input_strand}" != "no" ]; then
            >&2 echo "strand must be empty, 'reverse', 'yes', or 'no'"
            exit 1
        fi
    }

    output {
        File bam_after_input_validated = input_bam
    }
}