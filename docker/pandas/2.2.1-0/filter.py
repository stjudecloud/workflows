import csv
import pandas as pd
import numpy as np
import sys
import argparse

def get_args():
    parser = argparse.ArgumentParser(
        description="Filter probes based on standard deviation.")
    parser.add_argument(
        "--num_probes", type=int, default=10000, help="Number of probes to retain after filtering.")
    parser.add_argument(
        "beta", type=str, help="Beta values CSV file.")

    args = parser.parse_args()

    return args

if __name__ == "__main__":
    args = get_args()

    # Read beta values and compute standard deviation
    data = []
    with open("beta.csv", "r") as f:
        reader = csv.reader(f)
        header = next(reader)
        for i, line in enumerate(reader):
            probe = line[0]
            sd = np.std([float(x) for x in line[1:]])
            data.append([probe, sd])

    sd_df = pd.DataFrame(data, columns=["probe", "sd"]).set_index("probe")

    # Filter probes based on standard deviation
    filtered_probes = sd_df.sort_values("sd", ascending=False).head(args.num_probes).index

    pd.Series(filtered_probes, index=filtered_probes).to_csv("filtered_probes.csv", index=False)

    # Filter beta values
    with open("beta.csv", "r") as f:
        reader = csv.reader(f)
        header = next(reader)
        header[0] = "probe"
        beta_data = [line for line in reader if line[0] in filtered_probes]

    bd = pd.DataFrame(beta_data, columns=header).set_index("probe")
    bd.to_csv("filtered_beta.csv")