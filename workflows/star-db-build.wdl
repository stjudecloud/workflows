import "../tools/star.wdl"

workflow build_db {
    File reference_fasta 
    File gencode_gtf 
    String stardb_dir_name

    call star.build_db {
        input:
            reference_fasta=reference_fasta,
            gencode_gtf=gencode_gtf,
            stardb_dir_name=stardb_dir_name,
            ncpu=4
    }

    output {
        File out_dir = build_db.out_dir
    }
}