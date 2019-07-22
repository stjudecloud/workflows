task count {
    File bam
    File gff
    String strand = "reverse"
 
    command {
        htseq-count -f bam -r pos -s ${strand} -m union -i gene_name --secondary-alignments ignore --supplementary-alignments ignore ${bam} ${gff}
    }
}
