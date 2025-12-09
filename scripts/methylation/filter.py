import csv
import pandas as pd
import numpy as np
import argparse


def get_args():
    parser = argparse.ArgumentParser(
        description="Filter probes based on standard deviation."
    )
    parser.add_argument(
        "--num-probes",
        type=int,
        default=10000,
        help="Number of probes to retain after filtering.",
    )
    parser.add_argument("--output-name", type=str, help="Name for output file.")
    parser.add_argument(
        "--filtered-probes", type=str, help="Name for filtered probes file."
    )
    parser.add_argument(
        "--pval-threshold",
        type=float,
        default=0.01,
        help="P-value threshold for filtering.",
    )
    parser.add_argument(
        "--pval-sample-fraction",
        type=float,
        default=0.5,
        help="Fraction of samples that must exceed the p-value threshold.",
    )
    parser.add_argument("--pval", type=str, help="P-values CSV file.")
    parser.add_argument(
        "--exclude", type=str, nargs='*', help="Files with probes to exclude."
    )
    parser.add_argument("beta", type=str, help="Beta values CSV file.")
    args = parser.parse_args()

    return args


if __name__ == "__main__":
    args = get_args()

    # Read probes to exclude
    probes_to_exclude = []
    if args.exclude is not None:
        for exclude_file in args.exclude:
            with open(exclude_file, "r") as ef:
                for line in ef:
                    probes_to_exclude.append(line.strip())
        print("Number of probes to exclude:", len(probes_to_exclude))

    # Read p-values and find probes with high p-value in too many samples
    high_pval_probes = []
    if args.pval is not None:
        with open(args.pval, "r") as f:
            reader = csv.reader(f)
            header = next(reader)
            allowable_sample_count = args.pval_sample_fraction * (len(header) - 1)
            for i, line in enumerate(reader):
                if i % 10000 == 0:
                    print(f"Processing p-value for probe {i}")
                probe = line[0]
                # Get count of samples where the p-value exceeds the threshold
                count_high_pval = sum(
                    1 for x in line[1:] if float(x) > args.pval_threshold
                )
                # Check if the count exceeds the allowable fraction
                if count_high_pval > allowable_sample_count:
                    high_pval_probes.append(probe)
        print(
            "Number of probes with high p-value in too many samples:",
            len(high_pval_probes),
        )
        pd.Series(high_pval_probes).to_csv(
            "high_pval_probes.csv", index=False, header=False
        )

    # Combine with probes to exclude
    exclude_probe_list = set(high_pval_probes).union(set(probes_to_exclude))
    print("Total number of probes to exclude:", len(exclude_probe_list))

    # Read beta values and compute standard deviation
    data = []

    with open(args.beta, "r") as f:
        reader = csv.reader(f)
        header = next(reader)
        for i, line in enumerate(reader):
            if i % 10000 == 0:
                print(f"Processing probe {i}")

            probe = line[0]
            sd = np.std([float(x) for x in line[1:]])
            if probe not in exclude_probe_list:
                data.append([probe, sd])

    sd_df = pd.DataFrame(data, columns=["probe", "sd"]).set_index("probe")

    # Filter probes based on standard deviation
    filtered_probes = (
        sd_df.sort_values("sd", ascending=False).head(args.num_probes).index
    )

    pd.Series(filtered_probes, index=filtered_probes).to_csv(
        args.filtered_probes, index=False
    )

    # Filter beta values
    # NOTE: we are iterating over the file a second time to avoid loading the entire file into memory
    print("Reading beta for filtering")
    with open(args.beta, "r") as f:
        reader = csv.reader(f)
        header = next(reader)
        header[0] = "probe"
        beta_data = [line for line in reader if line[0] in filtered_probes]

    print("Writing filtered beta values")
    bd = pd.DataFrame(beta_data, columns=header).set_index("probe")
    bd.to_csv(args.output_name)
