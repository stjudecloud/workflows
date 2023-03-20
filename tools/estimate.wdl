## # ESTIMATE
##
## This WDL tool wraps the [ESTIMATE software package](https://bioinformatics.mdanderson.org/estimate/).

version 1.0

task calc_tpm {
    input {
        File counts
        File gene_lengths
        String outfile_name = basename(counts, ".feature-counts.txt") + ".TPM.txt"
        Int max_retries = 1
    }

    command <<<
        COUNTS="~{counts}" GENE_LENGTHS="~{gene_lengths}" OUTFILE="~{outfile_name}" python3 - <<END
import os  # lint-check: ignore

counts_file = open(os.environ['COUNTS'], 'r')
counts = {}
for line in counts_file:
    gene, count = line.split('\t')
    if gene[0:2] == '__':
        break
    counts[gene.strip()] = int(count.strip())
counts_file.close()

lengths_file = open(os.environ['GENE_LENGTHS'], 'r')
rpks = {}
tot_rpk = 0
lengths_file.readline()  # discard header
for line in lengths_file:
    gene, length = line.split('\t')
    rpk = counts[gene.strip()] / int(length.strip()) * 1000
    tot_rpk += rpk
    rpks[gene.strip()] = rpk
lengths_file.close()

sf = tot_rpk / 1000000

sample_name = '.'.join(os.environ['OUTFILE'].split('.')[:-2])
outfile = open(os.environ['OUTFILE'], 'w')
print(f"Gene name\t{sample_name}", file=outfile)
for gene, rpk in sorted(rpks.items()):
    tpm = rpk / sf
    print(f"{gene}\t{tpm:.3f}", file=outfile)
outfile.close()
END
    >>>

    runtime {
        memory: "4 GB"
        disk: "4 GB"
        docker: 'ghcr.io/stjudecloud/util:1.2.0'
        maxRetries: max_retries
    }

    output {
        File tpm_file = "~{outfile_name}"
    }
}

task run_ESTIMATE {
    input {
        File gene_expression_file
        String outfile_name = basename(gene_expression_file, ".TPM.txt") + ".ESTIMATE.gct"
        Int max_retries = 1
    }

    command <<<
        cp "~{gene_expression_file}" gene_expression.txt
        Rscript - <<END
library("estimate")

infile <- read.table(file = "gene_expression.txt", sep = '\t', header = TRUE)
filtered <- infile[infile$"Gene.name" %in% common_genes[['GeneSymbol']], ]
write.table(filtered, sep = "\t", file = "filtered.tsv", row.names = FALSE, quote = FALSE)
outputGCT("filtered.tsv", "gene_expression.gct")
estimateScore("gene_expression.gct", "common_estimate.gct", platform = "illumina")
END
    mv common_estimate.gct "~{outfile_name}"
    >>>

    runtime {
        memory: "4 GB"
        disk: "4 GB"
        docker: 'ghcr.io/stjudecloud/estimate:1.0.0'
        maxRetries: max_retries
    }

    output {
        File estimate_file = "~{outfile_name}"
    }
}