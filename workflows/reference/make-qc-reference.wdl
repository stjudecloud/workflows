version 1.0

import "../../tools/fastq_screen.wdl"

workflow make_qc_reference {
    call fastq_screen.build_db as fastq_screen_build_db

    output {
        File fastq_screen_db = fastq_screen_build_db.db
    }
}