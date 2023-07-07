## # ChIP-Seq Standard
##
## This WDL workflow runs the BWA ChIP-seq alignment workflow for St. Jude Cloud.
##
## The workflow takes an input BAM file and splits it into FastQ files for each read in the pair.
## The read pairs are then passed through BWA alignment to generate a BAM file.
## File validation is performed at several steps, including immediately preceeding output.
##
## ## LICENSING
##
## #### MIT License
##
## Copyright 2021-Present St. Jude Children's Research Hospital
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

version 1.1

# TODO should we enfore an order to imports? Alphabetical? Order referenced? URLs before/after relatives?
import "../general/bam-to-fastqs.wdl" as b2fq
import "https://raw.githubusercontent.com/stjude/seaseq/2.3/workflows/workflows/mapping.wdl" as seaseq_map
import "https://raw.githubusercontent.com/stjude/seaseq/3.0/workflows/tasks/seaseq_util.wdl" as seaseq_util
import "../../tools/ngsderive.wdl"
import "../../tools/picard.wdl"
import "../../tools/samtools.wdl"
import "https://raw.githubusercontent.com/stjude/seaseq/3.0/workflows/tasks/samtools.wdl" as seaseq_samtools
import "../../tools/util.wdl"
import "../../tools/deeptools.wdl"
import "../../tools/md5sum.wdl"

workflow chipseq_standard {
    parameter_meta {
        bam: "Input BAM format file to realign with bowtie"
        bowtie_indexes: "Database of v1 reference files for the bowtie aligner. Can be generated with https://github.com/stjude/seaseq/blob/master/workflows/tasks/bowtie.wdl. [*.ebwt]"
        excludelist: "Optional list of regions that will be excluded after reference alignment"
        prefix: "Prefix for output files"
        validate_input: "Run Picard ValidateSamFile on the input BAM"
        use_all_cores: "Use all cores for multi-core steps?"
        subsample_n_reads: "Only process a random sampling of `n` reads. <=`0` for processing entire input BAM."
        max_retries: "Number of times to retry failed steps"
    }

    input {
        File bam
        Array[File] bowtie_indexes
        File? excludelist
        String prefix = basename(bam, ".bam")
        Boolean validate_input = true
        Boolean use_all_cores = false
        Int subsample_n_reads = -1
        Int max_retries = 1
    }

    if (validate_input) {
        call picard.validate_bam as validate_input_bam { input:
            bam=bam,
            max_retries=max_retries
        }
    }

    if (subsample_n_reads > 0) {
        call samtools.subsample { input:
            bam=bam,
            desired_reads=subsample_n_reads,
            use_all_cores=use_all_cores,
            max_retries=max_retries
        }
    }
    File selected_bam = select_first([subsample.sampled_bam, bam])

    call util.get_read_groups { input:
        bam=selected_bam,
        format_for_star=false,
        max_retries=max_retries
    }
    Array[String] read_groups = read_lines(get_read_groups.read_groups_file)

    call b2fq.bam_to_fastqs { input:
        bam=selected_bam,
        paired_end=false,
        use_all_cores=use_all_cores,
        max_retries=max_retries
    }

    call samtools.index as samtools_index_input { input:
        bam=selected_bam,
        max_retries=max_retries
    }

    call ngsderive.read_length { input:
        bam=selected_bam,
        bam_index=samtools_index_input.bam_index,
        max_retries=max_retries
    }

    scatter (pair in zip(bam_to_fastqs.read1s, read_groups)){
        call seaseq_util.basicfastqstats as basic_stats { input:
            fastqfile=pair.left
        }
        call seaseq_map.mapping as bowtie_single_end_mapping { input:
            fastqfile=pair.left,
            index_files=bowtie_indexes,
            metricsfile=basic_stats.metrics_out,
            blacklist=excludelist
        }
        File chosen_bam = select_first(
            [
                bowtie_single_end_mapping.bklist_bam,
                bowtie_single_end_mapping.mkdup_bam,
                bowtie_single_end_mapping.sorted_bam
            ]
        )
        call util.add_to_bam_header { input:
            bam=chosen_bam,
            additional_header=pair.right,
            max_retries=max_retries
        }
        String rg_id_field = sub(sub(pair.right, ".*ID:", "ID:"), "\t.*", "") 
        String rg_id = sub(rg_id_field, "ID:", "")
        call samtools.addreplacerg as single_end { input:
            bam=add_to_bam_header.reheadered_bam,
            read_group_id=rg_id,
            max_retries=max_retries
        }
    }

    Array[File] aligned_bams = select_first([single_end.tagged_bam, []])
    scatter(aligned_bam in aligned_bams){
        call picard.clean_sam as picard_clean { input:
            bam=aligned_bam,
            max_retries=max_retries
        }
    }

    call picard.merge_sam_files as picard_merge { input:
        bams=picard_clean.cleaned_bam,
        prefix=prefix,
        max_retries=max_retries
    }

    call seaseq_samtools.markdup { input:
        bamfile=picard_merge.merged_bam,
        outputfile=prefix + ".bam"
    }
    call samtools.index as samtools_index { input:
        bam=markdup.mkdupbam,
        use_all_cores=use_all_cores,
        max_retries=max_retries
    }
    call picard.validate_bam { input: bam=markdup.mkdupbam, max_retries=max_retries }

    call md5sum.compute_checksum { input:
        file=markdup.mkdupbam,
        max_retries=max_retries
    }

    call deeptools.bam_coverage as deeptools_bam_coverage { input:
        bam=markdup.mkdupbam,
        bam_index=samtools_index.bam_index,
        prefix=prefix,
        max_retries=max_retries
    }

    output {
        File harmonized_bam = markdup.mkdupbam
        File bam_checksum = compute_checksum.md5sum
        File bam_index = samtools_index.bam_index
        File bigwig = deeptools_bam_coverage.bigwig
    }
}
