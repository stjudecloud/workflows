import "../tools/star.wdl"

workflow star_alignment {
    File read_one_fastq
    File read_two_fastq
    File stardb

    call star.alignment {
        input:
            read_one_fastq=read_one_fastq,
            read_two_fastq=read_two_fastq,
            stardb=stardb,
            outprefix="Output."
    }
}