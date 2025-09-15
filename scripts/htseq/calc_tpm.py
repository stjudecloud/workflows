def main(counts_name_path, gene_lengths_path, outfile_path, has_header):
    counts_file = open(counts_name_path, "r")
    counts = {}
    if has_header:
        counts_file.readline()
    for line in counts_file:
        gene, count = line.split("\t")
        if gene[0:2] == "__":
            break
        counts[gene.strip()] = int(count.strip())
    counts_file.close()

    lengths_file = open(gene_lengths_path, "r")
    rpks = {}  # Reads Per Kilobase
    tot_rpk = 0
    lengths_file.readline()  # discard header
    for line in lengths_file:
        gene, length = line.split("\t")
        rpk = counts[gene.strip()] / int(length.strip()) * 1000
        tot_rpk += rpk
        rpks[gene.strip()] = rpk
    lengths_file.close()

    scaling_factor = tot_rpk / 1000000

    sample_name = ".".join(outfile_path.split(".")[:-2])  # assumed to end in `.TPM.txt`
    outfile = open(outfile_path, "w")
    print(f"feature\t{sample_name}", file=outfile)
    for gene, rpk in sorted(rpks.items()):
        tpm = rpk / scaling_factor
        print(f"{gene}\t{tpm:.3f}", file=outfile)
    outfile.close()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("counts_path")
    parser.add_argument("gene_lengths_path")
    parser.add_argument("outfile_path")
    parser.add_argument("--counts_has_header", action="store_true")

    args = parser.parse_args()
    main(
        args.counts_path,
        args.gene_lengths_path,
        args.outfile_path,
        args.counts_has_header,
    )
