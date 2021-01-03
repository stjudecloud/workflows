version 1.0

task calc_tpm {
    input {
        File counts
        File gene_lengths
        String outfile = basename(counts, ".feature-counts.txt") + ".tpm.txt"
        Int max_retries = 1
    }

    command <<<
        COUNTS="~{counts}" GENE_LENGTHS="~{gene_lengths}" OUTFILE="~{outfile}" python3 - <<END
import os

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

outfile = open(os.environ['OUTFILE'], 'w')
print("Gene name\tTPM", file=outfile)
for gene, rpk in sorted(rpks.items()):
    tpm = rpk / sf
    print(f"{gene}\t{tpm:.3f}", file=outfile)
outfile.close()
END
    >>>

    runtime {
        memory: "4 GB"
        disk: "4 GB"
        docker: 'stjudecloud/util:branch-ESTIMATE-1.1.0'
        maxRetries: max_retries
    }

    output {
        File out = "~{outfile}"
    }
}