## **WARNING:** this workflow is experimental! Use at your own risk!
#
# SPDX-License-Identifier: MIT
# Copyright St. Jude Children's Research Hospital
version 1.1

import "../../tools/samtools.wdl"

workflow samtools_merge {
    input {
        Array[File] bams
        Int max_length = 100
        Int? max_retries
        String prefix = basename(bams[0], ".bam")
        Boolean use_all_cores = false
    }

    Int bam_length = length(bams)

    if (bam_length > max_length){
        # Find the number of merges required
        scatter ( merge_num in range((bam_length / max_length) + 1)){
            # Get the sublist of bams
            scatter ( bam_num in range(max_length)){
                Int num = if merge_num > 0 then bam_num + (merge_num * max_length) else bam_num
                if (num < bam_length){
                    File bam_list = bams[num]
                }
            }
        }
        scatter (list in bam_list){
            call samtools.merge as inner_merge { input:
                bams=select_all(list),
                prefix=prefix,
                combine_pg=false,
                use_all_cores=use_all_cores,
                max_retries=max_retries
            }        
        }
        call samtools.merge as final_merge { input:
            bams=inner_merge.merged_bam,
            prefix=prefix,
            attach_rg=false,
            combine_pg=true,
            combine_rg=true,
            use_all_cores=use_all_cores,
            max_retries=max_retries
        }
    }

    if (bam_length < max_length){
        call samtools.merge as basic_merge { input:
            bams=bams,
            prefix=prefix,
            combine_pg=false,
            use_all_cores=use_all_cores,
            max_retries=max_retries
        }
    }

    output {
        File merged_bam = select_first([final_merge.merged_bam, basic_merge.merged_bam])
    }
}
