version 1.0

import "../../tools/deeptools.wdl"
import "../../tools/samtools.wdl"
import "../../tools/picard.wdl"
import "../../tools/md5sum.wdl"
import "https://raw.githubusercontent.com/stjude/XenoCP/3.1.4/wdl/workflows/xenocp.wdl" as xenocp_workflow

workflow alignment_post {
    input {
        File bam
        File? contaminant_db
        String xenocp_aligner = ""
        Boolean cleanse_xenograft = false
        Boolean detect_nproc = false
        Int max_retries = 1
    }
    
    parameter_meta {
        detect_nproc: "Use all available cores for multi-core steps"
        max_retries: "Number of times to retry failed steps"
    }
    
    call picard.sort as picard_sort { input: bam=bam, max_retries=max_retries }
    
    call samtools.index as samtools_index { input:
        bam=picard_sort.sorted_bam,
        max_retries=max_retries,
        detect_nproc=detect_nproc
    }

    if (cleanse_xenograft){
        File contam_db = select_first([contaminant_db, ""])
        call xenocp_workflow.xenocp {
            input:
                input_bam=picard_sort.sorted_bam,
                input_bai=samtools_index.bai,
                reference_tar_gz=contam_db,
                aligner=xenocp_aligner,
                skip_duplicate_marking=true
        }
    }
    File aligned_bam = select_first([xenocp.bam, picard_sort.sorted_bam])
    File aligned_bai = select_first([xenocp.bam_index, samtools_index.bai])

    call picard.validate_bam { input: bam=aligned_bam, max_retries=max_retries }

    call md5sum.compute_checksum { input: infile=aligned_bam, max_retries=max_retries }
    
    call deeptools.bamCoverage as deeptools_bamCoverage { input:
        bam=aligned_bam,
        bai=aligned_bai,
        max_retries=max_retries
    }

    output {
        File out_bam = aligned_bam
        File bam_index = aligned_bai
        File bam_checksum = compute_checksum.outfile
        File bigwig = deeptools_bamCoverage.bigwig
    }
}