## scRNA-Seq Standard
##
## This WDL workflow runs the Cell Ranger scRNA-Seq alignment workflow for St. Jude Cloud.
##
## The workflow takes an input BAM file and splits it into FastQ files for each read in the pair. 
## The read pairs are then passed through Cell Ranger to generate a BAM file and perform
## quantification. Strandedness is inferred using ngsderive.
## File validation is performed at several steps, including immediately preceeding output.
##
## ## LICENSING
##
## #### MIT License
##
## Copyright 2022-Present St. Jude Children's Research Hospital
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

import "./10x-bam-to-fastqs.wdl" as b2fq
import "../../tools/picard.wdl"
import "../../tools/ngsderive.wdl"
import "../../tools/samtools.wdl"
import "../../tools/util.wdl"
import "../../tools/cellranger.wdl"
import "../../tools/md5sum.wdl"

workflow scrnaseq_standard {
    input {
        File bam
        String output_prefix = basename(bam, ".bam")
        File gtf
        File transcriptome_tar_gz
        String strandedness = ""
        Int subsample_n_reads = -1
        Boolean validate_input = true
        Boolean use_all_cores = false
        Int? max_retries
    }

    parameter_meta {
        bam: "Input BAM format file to quality check"
        output_prefix: "Prefix for output files"
        gtf: "Gzipped GTF feature file"
        transcriptome_tar_gz: "Database of reference files for Cell Ranger. Can be downloaded from 10x Genomics."
        strandedness: {
            description: "Strandedness protocol of the RNA-Seq experiment. If unspecified, strandedness will be inferred by `ngsderive`."
            choices: [
                '',
                'Stranded-Reverse',
                'Stranded-Forward',
                'Unstranded'
            ]
        },
        subsample_n_reads: "Only process a random sampling of `n` reads. <=`0` for processing entire input BAM."
        use_all_cores: "Use all cores for multi-core steps?"
        max_retries: "Number of times to retry failed steps. Overrides task level defaults."
    }

    String provided_strandedness = strandedness

    call parse_input { input: input_strand=provided_strandedness }
    if (validate_input) {
       call picard.validate_bam as validate_input_bam { input: bam=bam, max_retries=max_retries }
    }

    if (subsample_n_reads > 0) {
        call samtools.subsample {
            input:
                bam=bam,
                max_retries=max_retries,
                desired_reads=subsample_n_reads,
                use_all_cores=use_all_cores
        }
    }
    File selected_input_bam = select_first([subsample.sampled_bam, bam])

    call b2fq.cell_ranger_bam_to_fastqs { input: bam=selected_input_bam, max_retries=max_retries, use_all_cores=use_all_cores }

    call cellranger.count {
        input:
            fastqs_tar_gz=cell_ranger_bam_to_fastqs.fastqs_archive,
            transcriptome_tar_gz=transcriptome_tar_gz,
            id=output_prefix,
            max_retries=max_retries,
            use_all_cores=use_all_cores
    }
    call picard.validate_bam { input: bam=count.bam, max_retries=max_retries }
    call ngsderive.infer_strandedness as ngsderive_strandedness { input: bam=count.bam, bam_index=count.bam_index, gtf=gtf, max_retries=max_retries }

    call md5sum.compute_checksum { input: infile=count.bam, max_retries=max_retries }

    output {
        File harmonized_bam = count.bam
        File bam_checksum = compute_checksum.md5sum
        File bam_index = count.bam_index
        File qc = count.qc
        File barcodes = count.barcodes
        File features = count.features
        File matrix = count.matrix
        File filtered_gene_h5 = count.filtered_gene_h5
        File raw_gene_h5 = count.raw_gene_h5
        File raw_barcodes = count.raw_barcodes
        File raw_features = count.raw_features
        File raw_matrix = count.raw_matrix
        File mol_info_h5 = count.mol_info_h5
        File web_summary = count.web_summary
        File inferred_strandedness = ngsderive_strandedness.strandedness_file
    }
}

task parse_input {
    input {
        String input_strand
    }

    command {
        if [ -n "~{input_strand}" ] && [ "~{input_strand}" != "Stranded-Reverse" ] && [ "~{input_strand}" != "Stranded-Forward" ] && [ "~{input_strand}" != "Unstranded" ]; then
            >&2 echo "strandedness must be empty, 'Stranded-Reverse', 'Stranded-Forward', or 'Unstranded'"
            exit 1
        fi
    }

    runtime {
        memory: "4 GB"
        disk: "1 GB"
        docker: 'ghcr.io/stjudecloud/util:1.2.0'
    }

    output {
        String input_check = "passed"
    }
}
