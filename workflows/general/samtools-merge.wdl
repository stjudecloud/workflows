## **WARNING:** this workflow is experimental! Use at your own risk!

version 1.1

import "../../tools/samtools.wdl"

workflow samtools_merge {
    meta{
        description: "Runs `samtools merge`, with optional iteration to avoid maximum command line argument length"
        category: "Utility"
        outputs: {
            merged_bam: "The BAM resulting from merging all the input BAMs"
        }
        allowNestedInputs: true
    }

    parameter_meta{
        bams: "BAMs to merge into a final BAM"
        prefix: "Prefix for output BAM."
        use_all_cores: "Use all cores? Recommended for cloud environments."
        max_length: "Maximum number of BAMs to merge before using iteration"
    }

    input {
        Array[File] bams
        String prefix
        Boolean use_all_cores = false
        Int max_length = 100
    }

    Int bam_length = length(bams)

    if (bam_length > max_length){
        # Find the number of merges required
        scatter (merge_num in range((bam_length / max_length) + 1)){
            # Get the sublist of bams
            scatter (bam_num in range(max_length)){
                Int num = (
                    if merge_num > 0
                    then bam_num + (merge_num * max_length)
                    else bam_num
                )
                if (num < bam_length){
                    File bam_list = bams[num]
                }
            }
        }
        scatter (list in bam_list){
            call samtools.merge as inner_merge { input:
                bams = select_all(list),
                prefix,
                attach_rg = false,
                combine_pg = false,
                use_all_cores,
            }
        }
        call samtools.merge as final_merge { input:
            bams = inner_merge.merged_bam,
            prefix,
            attach_rg = false,
            combine_pg = true,
            combine_rg = true,
            use_all_cores,
        }
    }

    if (bam_length < max_length){
        call samtools.merge as basic_merge { input:
            bams,
            prefix,
            attach_rg = false,
            combine_pg = false,
            use_all_cores,
        }
    }

    output {
        File merged_bam = select_first([final_merge.merged_bam, basic_merge.merged_bam])
    }
}
