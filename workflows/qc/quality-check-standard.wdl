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

import "https://raw.githubusercontent.com/stjudecloud/workflows/master/tools/samtools.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/master/tools/picard.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/master/tools/qc.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/master/tools/fastqc.wdl" as fqc
import "https://raw.githubusercontent.com/stjudecloud/workflows/master/tools/qualimap.wdl"
#import "https://raw.githubusercontent.com/stjudecloud/workflows/master/tools/ngsderive.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/master/tools/deeptools.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/master/tools/md5sum.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/master/tools/multiqc.wdl" as mqc

workflow quality_check {
    input {
        File bam
        Int max_retries = 1
    }

    call samtools.quickcheck { input: bam=bam, max_retries=max_retries }
    call picard.validate_bam { input: bam=bam, max_retries=max_retries }
    call qc.parse_validate_bam { input: in=validate_bam.out, max_retries=max_retries }
    call qualimap.bamqc as qualimap_bamqc { input: bam=bam, max_retries=max_retries }
#    call qualimap.rnaseq as qualimap_rnaseq { input: bam=bam, gencode_gtf=gencode_gtf, strand=strand, inferred=ngsderive_strandedness.strandedness, max_retries=max_retries }
    call fqc.fastqc { input: bam=bam, max_retries=max_retries }
    call samtools.flagstat as samtools_flagstat { input: bam=bam, max_retries=max_retries }
    call samtools.index as samtools_index { input: bam=bam, max_retries=max_retries }
    call deeptools.bamCoverage as deeptools_bamCoverage { input: bam=bam, bai=samtools_index.bai, max_retries=max_retries }
#    call ngsderive.infer_strand as ngsderive_strandedness { input: bam=bam, bai=samtools_index.bai, gtf=gencode_gtf, max_retries=max_retries }
    call md5sum.compute_checksum { input: infile=bam, max_retries=max_retries }
    call mqc.multiqc {
        input:
            sorted_bam=bam,
            validate_sam_string=validate_bam.out,
            qualimap_bamqc=qualimap_bamqc.out_files,
#            qualimap_rnaseq=qualimap_rnaseq.out_files,
            fastqc_files=fastqc.out_files,
            flagstat_file=samtools_flagstat.outfile,
            bigwig_file=deeptools_bamCoverage.bigwig,
            max_retries=max_retries
    }

    output {
        File bam_checksum = compute_checksum.outfile
        File bam_index = samtools_index.bai
        File bigwig = deeptools_bamCoverage.bigwig
        File flagstat = samtools_flagstat.outfile
        Array[File] fastqc_results = fastqc.out_files
        Array[File] qualimap_bamqc_results = qualimap_bamqc.out_files
#        Array[File] qualimap_rnaseq_results = qualimap_rnaseq.out_files
        File multiqc_zip = multiqc.out
#        File inferred_strandedness = ngsderive_strandedness.strandedness_file
    }
}
