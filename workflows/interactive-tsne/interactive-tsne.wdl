## Description:
##
## This WDL workflow remaps input bams and compares them to St. Jude RNA-seq samples using 
## a t-SNE plot. 
##
## Inputs:
##
## in_bams - bams to compare to reference data
## tissue_type - tissue type of the input samples [blood, brain, solid]
## gencode_gtf - Gencode annotation file
##
## Output:
##
## tsne_plot - t-SNE plot including the input samples compared to reference data
## generated_counts - RNAseq count data for the input bams
## generated_mappings - Bam output from the St. Jude Cloud RNA-seq workflow for each input bam
##
## LICENSING:
## 
## MIT License
##
## Copyright 2020 St. Jude Children's Research Hospital
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

import "https://raw.githubusercontent.com/stjudecloud/workflows/master/workflows/rnaseq/rnaseq-standard.wdl" as rnav2
import "https://raw.githubusercontent.com/stjudecloud/workflows/master/tools/util.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/master/tools/tsne.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/master/tools/tar.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/master/tools/gzip.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/master/tools/wget.wdl"

workflow interactive_tsne {
    input { 
        Array[File] in_bams
        File gencode_gtf
        File? stardb_tar_gz
        File? reference_counts
        File? covariates_file
        File? gene_blacklist
        String tissue_type
        String? output_filename
    }
    
    if(! defined(gene_blacklist)){
        String blacklist_url = "https://stjudecloud.blob.core.windows.net/interactive-tsne/gene.blacklist.tsv"
        call wget.download as blacklist_download { input: url=blacklist_url, outfilename="blacklist.txt" }
    }

    File blacklist = select_first([gene_blacklist, blacklist_download.outfile])

    String outfile = select_first([output_filename, "output.html"])
    
    if (! defined(stardb_tar_gz)){
        String star_url = "https://stjudecloud.blob.core.windows.net/reference/RNA-seq/v2/STARDB.tar.gz"
        call wget.download as star_download { input: url=star_url, outfilename="STARDB.tar.gz" }
    }
    File stardb = select_first([stardb_tar_gz, star_download.outfile])

    call tsne.validate_tissue_type { input: tissue_type=tissue_type }
 
    call gzip.unzip as uncompress_gencode { input: infile=gencode_gtf }
    scatter (bam in in_bams) {
        String name = basename(bam, ".bam")
        call rnav2.rnaseq_standard { input: gencode_gtf=uncompress_gencode.outfile, input_bam=bam, stardb_tar_gz=stardb, output_prefix=name + '.' }
    }

    if (! defined(reference_counts)){
        # Based on tissue type selection we can provide reference data
        String reference_url = if (tissue_type == 'blood') then "https://stjudecloud.blob.core.windows.net/interactive-tsne/blood_counts.tar.gz" else
                               if (tissue_type == 'brain') then "https://stjudecloud.blob.core.windows.net/interactive-tsne/brain_counts.tar.gz" else
                               if (tissue_type == 'solid') then "https://stjudecloud.blob.core.windows.net/interactive-tsne/solid_counts.tar.gz"
                               else "" 
        call wget.download as reference_counts_download { input: url=reference_url, outfilename="counts.tar.gz"}

        String covariates_url = if (tissue_type == 'blood') then "https://stjudecloud.blob.core.windows.net/interactive-tsne/blood_covariates.test.tsv" else
                               if (tissue_type == 'brain') then "https://stjudecloud.blob.core.windows.net/interactive-tsne/brain_covariates.test.tsv" else
                               if (tissue_type == 'solid') then "https://stjudecloud.blob.core.windows.net/interactive-tsne/solid_covariates.test.tsv"
                               else "" 
        call wget.download as covariates_download { input: url=covariates_url, outfilename="covariates.tsv" }
    }

    File reference_file = select_first([reference_counts, reference_counts_download.outfile])
    File covariates_input = select_first([covariates_file, covariates_download.outfile])

    call gzip.unzip { input: infile=reference_file }
    call tar.untar { input: infile=unzip.outfile }

    scatter (bam in in_bams){
        call util.file_prefix { input: in_file=bam }
    }

    # generate covariates information for input samples and add to covariates file.
    
    call tsne.append_input { input: inputs=file_prefix.out, covariates_infile=covariates_input}
    
    call tsne.plot as generate_plot{
        input:
            counts=untar.outfiles,
            inputs=file_prefix.out,
            input_counts=rnaseq_standard.gene_counts,
            blacklist=blacklist,
            covariates=append_input.covariates_outfile,
            gencode_gtf=gencode_gtf,
            outfile=outfile
    }
     
    output {
        File tsne_plot = generate_plot.html
        Array[File] generated_counts = rnaseq_standard.gene_counts
        Array[File] generated_mappings = rnaseq_standard.bam
    } 
}
