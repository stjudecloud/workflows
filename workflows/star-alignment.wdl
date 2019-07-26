import "../tools/star.wdl"

workflow star_alignment {
    Array[File] read_one_fastqs
    Array[File] read_two_fastqs
    File stardb_dir
    String output_prefix
    String? read_groups

    call star.alignment {
        input:
            read_one_fastqs=read_one_fastqs,
            read_two_fastqs=read_two_fastqs,
            stardb_dir=stardb_dir,
            output_prefix=output_prefix,
            read_groups=read_groups
    }

    output {
       File bam = alignment.bam 
    }
}
