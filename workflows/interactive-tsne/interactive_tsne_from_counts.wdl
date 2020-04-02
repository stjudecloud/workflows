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

import "https://raw.githubusercontent.com/stjudecloud/workflows/master/tools/util.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/master/tools/tsne.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/master/tools/tar.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/master/tools/gzip.wdl"
import "https://raw.githubusercontent.com/stjudecloud/workflows/master/tools/wget.wdl"

workflow interactive_tsne_from_counts {
    input { 
        Array[File] in_counts
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

    call gzip.unzip { input: infile=reference_file }
    call tar.untar { input: infile=unzip.outfile }

    scatter (count in in_counts){
        call util.file_prefix { input: in_file=count }
    }

    # generate covariates information for input samples and add to covariates file.   
    call tsne.append_input { input: inputs=file_prefix.out, covariates_infile=covariates_input}
    
    call tsne.plot as generate_plot{
        input:
            counts=untar.outfiles,
            inputs=file_prefix.out,
            input_counts=in_counts,
            blacklist=gene_blacklist,
            covariates=append_input.covariates_outfile,
            gencode_gtf=gencode_gtf,
            outfile=output_filename,
            tissue_type=tissue_type
    }
     
    output {
        File tsne_plot = generate_plot.html
    }

    parameter_meta {
        in_counts: {
            help: "Provide count files to run for comparison"
        }
        tissue_type: {
            help: "Provide the tissue type to compare against: [blood, brain, solid]",
            choices: ['blood', 'brain', 'solid']
        }
        output_filename: "Name for the output HTML t-SNE plot"
    }
}
