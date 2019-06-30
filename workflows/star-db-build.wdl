import "../tools/star.wdl"

workflow build_db {
    File reference_fasta 
    File gencode_gtf 
    File stardb

    call star.build_db {
        input:
            reference_fasta=reference_fasta,
            gencode_gtf=gencode_gtf,
            stardb=stardb,
            ncpu=4
    }
}