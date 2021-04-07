## # RNA-Seq Expression Classification from BAMs
##
## This WDL workflow compares gene counts from BAM files to St. Jude RNA-seq samples using 
## a t-SNE plot.
##
## * The BAM-based RNA-Seq Expression Classification pipeline requires that alignment must be against the [GRCh38_no_alt reference](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz). It should use parameters as specifed in our [RNA-seq workflow](https://stjudecloud.github.io/rfcs/0001-rnaseq-workflow-v2.0.0.html) to minimize any discrepancies caused by differing alignment specification.
##
## ### Output
##
## tsne_plot
## : t-SNE plot including the input samples compared to reference data
##
## generated_counts
## : RNAseq count data for the input bams
##
## ## LICENSING
## 
## #### MIT License
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

import "https://raw.githubusercontent.com/stjudecloud/workflows/rnaseq-standard/v3.0.0/tools/htseq.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/rnaseq-standard/v3.0.0/tools/samtools.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/rnaseq-standard/v3.0.0/tools/ngsderive.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/rnaseq-standard/v3.0.0/tools/util.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/rnaseq-standard/v3.0.0/tools/tsne.wdl"

workflow rnaseq_expression_classification_from_bams {
    input { 
        Array[File] in_bams
        File gencode_gtf
        File? reference_counts
        File? covariates_file
        File gene_blacklist
        String tissue_type
        String output_filename = "output.html"
        File? blood_counts
        File? brain_counts
        File? solid_counts
        File? blood_covariates
        File? brain_covariates
        File? solid_covariates
    }
 
    scatter (bam in in_bams) {       
        call samtools.index as index { input: bam=bam}
        call ngsderive.infer_strandedness as infer { input: bam=bam, bai=index.bai, gtf=gencode_gtf}
        call htseq.count as count { input: bam=bam, gtf=gencode_gtf, provided_strandedness="", inferred_strandedness=infer.strandedness}
    }

    if (! defined(reference_counts)){
        # Based on tissue type selection we can provide reference data
        File? reference = if (tissue_type == 'blood') then blood_counts else
                               if (tissue_type == 'brain') then brain_counts else
                               if (tissue_type == 'solid') then solid_counts
                               else "" 
        
        File? covariates = if (tissue_type == 'blood') then blood_covariates else
                               if (tissue_type == 'brain') then brain_covariates else
                               if (tissue_type == 'solid') then solid_covariates
                               else "" 
        
    }

    File reference_file = select_first([reference_counts, reference])
    File covariates_input = select_first([covariates_file, covariates])

    scatter (bam in in_bams){
        call util.file_prefix { input: in_file=bam }
    }

    # generate covariates information for input samples and add to covariates file.   
    call tsne.append_input { input: inputs=file_prefix.out, covariates_infile=covariates_input}
    
    call tsne.plot as generate_plot{
        input:
            counts=reference_file,
            inputs=file_prefix.out,
            input_counts=count.out,
            blacklist=gene_blacklist,
            covariates=append_input.covariates_outfile,
            gencode_gtf=gencode_gtf,
            outfile=output_filename,
            tissue_type=tissue_type
    }
     
    output {
        File tsne_plot = generate_plot.html
        File tsne_matrix = generate_plot.matrix
        Array[File] generated_counts = count.out
    }

    parameter_meta {
        in_bams: {
            help: "Provide bams to run for comparison",
            patterns: ["*.bam"]
        }
        tissue_type: {
            help: "Provide the tissue type to compare against",
            choices: ['blood', 'brain', 'solid']
        }
        output_filename: "Name for the output HTML RNA-Seq Expression Classification plot"
    }
}
