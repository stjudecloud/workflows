version 1.0

import "https://raw.githubusercontent.com/stjudecloud/workflows/build-db-refactor/tools/util.wdl"

workflow calc_gene_lengths {
    input {
        File gtf
    }

    call util.calc_gene_lengths as calc { input: gtf=gtf }

    output {
        File gene_lengths=calc.out
    }
}
