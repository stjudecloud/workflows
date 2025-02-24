import argparse


def get_args():
    parser = argparse.ArgumentParser(description="Filter Hi-C data.")
    parser.add_argument("--prefix", type=str, help="Prefix for output file.")
    parser.add_argument("--all_valid_pairs", type=str, help="All valid pairs file.")
    parser.add_argument("--filter_pairs", type=str, help="Filter pairs file.")

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = get_args()

    exclude_ID = set()
    with open(args.filter_pairs) as file:
        for line in file:
            exclude_ID.add(line.strip().split("\t")[0])

    outfile1 = args.prefix + ".allValidPairs.filtered"
    outfile2 = args.prefix + ".allValidPairs.removed"
    with open(args.all_valid_pairs, 'r') as pairs_file, open(outfile1, 'w') as of1, open(outfile2, 'w') as of2:
        for line in pairs_file:
            if line.strip().split("\t")[0] in exclude_ID:
                of2.write(line)
            else:
                of1.write(line)
