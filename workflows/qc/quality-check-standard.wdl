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

import "https://raw.githubusercontent.com/stjudecloud/workflows/qualimap-tar/tools/samtools.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/qualimap-tar/tools/picard.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/qualimap-tar/tools/qc.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/qualimap-tar/tools/qualimap.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/qualimap-tar/tools/ngsderive.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/qualimap-tar/tools/fastqc.wdl" as fqc
import "https://raw.githubusercontent.com/stjudecloud/workflows/qualimap-tar/tools/fastq_screen.wdl" as fq_screen
import "https://raw.githubusercontent.com/stjudecloud/workflows/qualimap-tar/tools/fq.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/qualimap-tar/tools/md5sum.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/qualimap-tar/tools/multiqc.wdl" as mqc

workflow quality_check {
    input {
        File gencode_gtf
        File fastq_screen_db
        String fastq_format = "sanger"
        String experiment
        File bam
        String strand = ""
        Int max_retries = 1
    }

    String provided_strand = strand

    call parse_input {
        input:
            input_experiment=experiment,
            input_strand=provided_strand,
            input_fq_db=fastq_screen_db,
            input_fq_format=fastq_format
    }

    call samtools.quickcheck { input: bam, max_retries=max_retries, wait_var=parse_input.input_check }
    call picard.validate_bam { input: bam, max_retries=max_retries, wait_var=parse_input.input_check }
    call qc.parse_validate_bam { input: in=validate_bam.out, max_retries=max_retries }
    call samtools.index as samtools_index { input: bam, max_retries=max_retries, wait_var=parse_input.input_check }
    call qualimap.bamqc as qualimap_bamqc { input: bam, max_retries=max_retries, wait_var=parse_input.input_check }
    call fqc.fastqc { input: bam, max_retries=max_retries, wait_var=parse_input.input_check }
    call samtools.flagstat as samtools_flagstat { input: bam, max_retries=max_retries, wait_var=parse_input.input_check }

    call samtools.subsample as samtools_subsample { input: bam, max_retries=max_retries, wait_var=parse_input.input_check }
    call picard.bam_to_fastq as b2fq { input: bam=samtools_subsample.sampled_bam, max_retries=max_retries }
    call fq.fqlint { input: read1=b2fq.read1, read2=b2fq.read2, max_retries=max_retries }
    call fq_screen.fastq_screen as fastq_screen { input: read1=b2fq.read1, read2=b2fq.read2, db=fastq_screen_db, format=fastq_format, max_retries=max_retries }
    
    call md5sum.compute_checksum { input: infile=bam, max_retries=max_retries }
    call ngsderive.instrument as ngsderive_instrument { input: bam, max_retries=max_retries, wait_var=parse_input.input_check }
    call ngsderive.readlen as ngsderive_readlen { input: bam, max_retries=max_retries, wait_var=parse_input.input_check }

    if (experiment == "RNA") {
        call ngsderive.infer_strand as ngsderive_strandedness { input: bam, bai=samtools_index.bai, gtf=gencode_gtf, max_retries=max_retries }
        call qualimap.rnaseq as qualimap_rnaseq { input: bam, gencode_gtf=gencode_gtf, provided_strand=provided_strand, inferred_strand=ngsderive_strandedness.strandedness, max_retries=max_retries }
        call mqc.multiqc as multiqc_rnaseq {
            input:
                sorted_bam=wait_var=parse_input.input_check,
                validate_sam_string=validate_bam.out,
                qualimap_bamqc=qualimap_bamqc.results,
                qualimap_rnaseq=qualimap_rnaseq.results,
                fastqc_files=fastqc.out_files,
                fastq_screen=fastq_screen.out_files,
                flagstat_file=samtools_flagstat.outfile,
                max_retries=max_retries
        }
    }
    if (experiment != "RNA") {
        call mqc.multiqc {
            input:
                sorted_bam=bam,
                validate_sam_string=validate_bam.out,
                qualimap_bamqc=qualimap_bamqc.results,
                fastqc_files=fastqc.out_files,
                fastq_screen=fastq_screen.out_files,
                flagstat_file=samtools_flagstat.outfile,
                max_retries=max_retries
        }
    }

    output {
        File bam_checksum = compute_checksum.outfile
        File bam_index = samtools_index.bai
        File flagstat = samtools_flagstat.outfile
        Array[File] fastqc_results = fastqc.out_files
        File qualimap_bamqc_results = qualimap_bamqc.results
        File? qualimap_rnaseq_results = qualimap_rnaseq.results
        Array[File] fastq_screen_results = fastq_screen.out_files
        File? multiqc_zip = multiqc.out
        File? multiqc_rnaseq_zip = multiqc_rnaseq.out
        File? inferred_strandedness = ngsderive_strandedness.strandedness_file
        File instrument_file = ngsderive_instrument.instrument_file
        File readlen_file = ngsderive_readlen.readlen_file
    }
}

task parse_input {
    input {
        String input_experiment
        String input_strand
        File input_fq_db
        String input_fq_format
    }

    command {
        if [ "~{input_experiment}" != "WGS" ] && [ "~{input_experiment}" != "WES" ] && [ "~{input_experiment}" != "RNA" ]; then
            >&2 echo "experiment input must be 'WGS', 'WES', or 'RNA'"
            exit 1
        fi

        if [ -n "~{input_strand}" ] && [ "~{input_strand}" != "stranded-reverse" ] && [ "~{input_strand}" != "stranded-forward" ] && [ "~{input_strand}" != "unstranded" ]; then
            >&2 echo "strand must be empty, 'stranded-reverse', 'stranded-forward', or 'unstranded'"
            exit 1
        fi

        if [ "$(basename ~{input_fq_db})" != "fastq-screen-db.tar.gz" ]; then
            >&2 echo "FastQ Screen database (input \"fastq_screen_db\") must be archived and named fastq-screen-db.tar.gz"
            exit 1
        fi

        if [ -n "~{input_fq_format}" ] && [ "~{input_fq_format}" != "sanger" ] && [ "~{input_fq_format}" != "illunima1.3" ]; then
            >&2 echo "fastq_format must be empty, 'sanger', or 'illumina1.3'"
            exit 1
        fi
    }

    output {
        String input_check = ""
    }
}