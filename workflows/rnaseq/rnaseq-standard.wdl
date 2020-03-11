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

import "https://raw.githubusercontent.com/stjudecloud/workflows/rfcs/qc-workflow/workflows/general/bam-to-fastqs.wdl" as b2fq
import "https://raw.githubusercontent.com/stjudecloud/workflows/rfcs/qc-workflow/tools/star.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/rfcs/qc-workflow/tools/picard.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/rfcs/qc-workflow/tools/fastqc.wdl" as fqc
import "https://raw.githubusercontent.com/stjudecloud/workflows/rfcs/qc-workflow/tools/ngsderive.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/rfcs/qc-workflow/tools/qualimap.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/rfcs/qc-workflow/tools/htseq.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/rfcs/qc-workflow/tools/samtools.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/rfcs/qc-workflow/tools/md5sum.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/rfcs/qc-workflow/tools/multiqc.wdl" as mqc
import "https://raw.githubusercontent.com/stjudecloud/workflows/rfcs/qc-workflow/tools/qc.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/rfcs/qc-workflow/tools/util.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/rfcs/qc-workflow/tools/deeptools.wdl"

workflow rnaseq_standard {
    input {
        File gencode_gtf
        File input_bam
        File stardb_tar_gz
        String strand = ""
        String output_prefix = "out"
        Int max_retries = 1
    }

    String provided_strand = strand

    call parse_input { input: input_strand=provided_strand }

    call samtools.get_read_groups { input: bam=input_bam, max_retries=max_retries, wait_var=parse_input.input_check }
    call util.prepare_read_groups_for_star { input: read_groups=get_read_groups.out, max_retries=max_retries }
    call b2fq.bam_to_fastqs { input: bam=input_bam, max_retries=max_retries, wait_var=parse_input.input_check }
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
    call ngsderive.infer_strand as ngsderive_strandedness { input: bam=picard_sort.sorted_bam, bai=samtools_index.bai, gtf=gencode_gtf, max_retries=max_retries }
    call htseq.count as htseq_count { input: bam=picard_sort.sorted_bam, gtf=gencode_gtf, provided_strand=provided_strand, inferred_strand=ngsderive_strandedness.strandedness, max_retries=max_retries }
    call md5sum.compute_checksum { input: infile=picard_sort.sorted_bam, max_retries=max_retries }

    output {
        File bam = picard_sort.sorted_bam
        File bam_checksum = compute_checksum.outfile
        File bam_index = samtools_index.bai
        File gene_counts = htseq_count.out
        File inferred_strandedness = ngsderive_strandedness.strandedness_file
    }
}

task parse_input {
    input {
        String input_strand
    }

    command {
        if [ -n "~{input_strand}" ] && [ "~{input_strand}" != "reverse" ] && [ "~{input_strand}" != "yes" ] && [ "~{input_strand}" != "no" ]; then
            >&2 echo "strand must be empty, 'reverse', 'yes', or 'no'"
            exit 1
        fi
    }

    output {
        String input_check = "passed"
    }
}