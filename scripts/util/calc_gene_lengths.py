import gtfparse
import numpy as np
from collections import defaultdict


def main(gtf_path, outfile_path, id_attr):
    gtf = gtfparse.read_gtf(gtf_path)

    only_exons = gtf[gtf["feature"] == "exon"]
    exon_starts = defaultdict(lambda: [])
    exon_ends = defaultdict(lambda: [])
    gene_start_offset = {}
    gene_end_offset = {}
    gene_exon_intersection = {}

    for (
        _index,
        value,
    ) in only_exons.iterrows():  # pyright: ignore [reportAttributeAccessIssue]
        feature_id = value[id_attr]
        start = value["start"]
        end = value["end"] + 1  # end is inclusive in GTF
        exon_starts[feature_id].append(start)
        exon_ends[feature_id].append(end)
        if feature_id not in gene_start_offset:
            gene_start_offset[feature_id] = start
            gene_end_offset[feature_id] = end
        else:
            gene_start_offset[feature_id] = min(
                gene_start_offset[feature_id], start
            )  # pyright: ignore [reportArgumentType]
            gene_end_offset[feature_id] = max(gene_end_offset[feature_id], end)

    for feature_id in exon_starts:
        gene_exon_intersection[feature_id] = np.full(
            gene_end_offset[feature_id] - gene_start_offset[feature_id], False
        )

        for start, end in zip(exon_starts[feature_id], exon_ends[feature_id]):
            gene_exon_intersection[feature_id][
                start
                - gene_start_offset[feature_id] : end
                - gene_start_offset[feature_id]
            ] = True

    outfile = open(outfile_path, "w")
    print("feature\tlength", file=outfile)
    for gene, exonic_intersection in sorted(gene_exon_intersection.items()):
        # np.count_nonzero() is faster than sum
        # np.count_nonzero() evaluates the "truthfulness" of
        # of all elements (by calling their '.__bool__()' method)
        length = np.count_nonzero(exonic_intersection)
        print(f"{gene}\t{length}", file=outfile)

    outfile.close()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("gtf_path", type=str)
    parser.add_argument("outfile_path", type=str)
    parser.add_argument("--id_attr", type=str, default="gene_name")

    args = parser.parse_args()
    main(args.gtf_path, args.outfile_path, args.id_attr)
