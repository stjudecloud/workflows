import sys
import argparse
import pandas as pd


def get_args():
    parser = argparse.ArgumentParser(description="Combine CSV files.")
    parser.add_argument("--output-name", type=str, help="Name for output file.")
    parser.add_argument("csvs", type=str, nargs="+", help="List of CSV files.")

    args = parser.parse_args()

    return args


def read(filename):
    df = pd.read_csv(filename)
    df.rename(
        columns={df.columns[0]: "probe"}, inplace=True
    )  # pyright: ignore [reportUnhashable]
    df.set_index("probe", inplace=True)
    return df


if __name__ == "__main__":
    args = get_args()

    # Combine data
    df = pd.concat([read(f) for f in args.csvs], axis=1, join="inner")
    df.to_csv(args.output_name, index=True)
