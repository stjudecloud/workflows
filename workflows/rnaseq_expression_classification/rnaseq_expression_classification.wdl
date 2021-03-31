## # RNA-Seq Expression Classification
##
## This WDL workflow remaps input bams and compares them to St. Jude RNA-seq samples using 
## a t-SNE plot. 
##
## ### Output
##
## tsne_plot
## : t-SNE plot including the input samples compared to reference data
##
## generated_counts
## : RNAseq count data for the input bams
##
## generated_mappings
## : Bam output from the St. Jude Cloud RNA-seq workflow for each input bam
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

import "https://raw.githubusercontent.com/stjudecloud/workflows/rnaseq-standard/v2.0.0/workflows/rnaseq/rnaseq-standard.wdl" as rnav2
import "https://raw.githubusercontent.com/stjudecloud/workflows/rnaseq-standard/v2.0.0/tools/util.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/rnaseq-standard/v2.0.0/tools/tsne.wdl"

workflow rnaseq_expression_classification {
    input { 
        Array[File] in_bams
        File gencode_gtf
        File stardb_tar_gz
        File? reference_counts
        File? covariates_file
        File gene_blacklist
        String tissue_type
        String output_filename = "output.html"
        String strandedness = ""
        Int subsample_n_reads = -1
        File? blood_counts
        File? brain_counts
        File? solid_counts
        File? blood_covariates
        File? brain_covariates
        File? solid_covariates
        Boolean detect_nproc = false
    }

    scatter (bam in in_bams) {
        String name = basename(bam, ".bam")
        call rnav2.rnaseq_standard {
            input:
                gencode_gtf=gencode_gtf,
                input_bam=bam,
                stardb_tar_gz=stardb_tar_gz,
                output_prefix=name,
                detect_nproc=detect_nproc,
                strandedness=strandedness, 
                subsample_n_reads=subsample_n_reads
        }
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
            input_counts=rnaseq_standard.gene_counts,
            blacklist=gene_blacklist,
            covariates=append_input.covariates_outfile,
            gencode_gtf=gencode_gtf,
            outfile=output_filename,
            tissue_type=tissue_type
    }
     
    output {
        File tsne_plot = generate_plot.html
        File tsne_matrix = generate_plot.matrix
        Array[File] generated_counts = rnaseq_standard.gene_counts
        Array[File] generated_mappings = rnaseq_standard.bam
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
        detect_nproc: "Determine number of cores and use all for multi-core steps"
        strandedness: {
            help: "Choose strandedness for samples. If omitted, the strandedness will attempt to be inferred.",
            choices: ['','Stranded-Forward', 'Stranded-Reverse', 'Unstranded'],
            default: ''
        }
    }
}
