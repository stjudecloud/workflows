import argparse
import pandas as pd


def get_args():
    parser = argparse.ArgumentParser(description="Combine CSV files.")
    parser.add_argument(
        "--output-name", type=str, help="Name for output file.", default="combined.csv"
    )
    parser.add_argument(
        "--simple-merge",
        action="store_true",
        help="Use simple merge rather than batched read. Use this if different arrays are to be combined.",
    )
    parser.add_argument("csvs", type=str, nargs="+", help="List of CSV files.")

    args = parser.parse_args()

    return args


def read(filename):
    df = pd.read_csv(filename)
    df.rename(
        columns={df.columns[0]: "probe"},  # pyright: ignore [reportUnhashable]
        inplace=True,
    )
    df.set_index("probe", inplace=True)
    return df


if __name__ == "__main__":
    args = get_args()

    # Combine data
    if args.simple_merge:
        dfs_to_concat = []
        for f in args.csvs:
            df = read(f)
            dfs_to_concat.append(df)
        df = pd.concat(dfs_to_concat, axis=1, join="inner")
        df.to_csv(args.output_name, index=True)
    else:
        dfs = []
        for f in args.csvs:
            dfs.append(pd.read_csv(f, chunksize=50000))

        for index, chunk in enumerate(dfs[0]):
            dfs_to_concat = []
            chunk = chunk.rename(columns={chunk.columns[0]: "probe"})
            chunk.set_index("probe", inplace=True)
            dfs_to_concat.append(chunk)
            for df in dfs[1:]:
                next_chunk = df.get_chunk()
                next_chunk = next_chunk.rename(columns={next_chunk.columns[0]: "probe"})
                next_chunk.set_index("probe", inplace=True)
                dfs_to_concat.append(next_chunk)
            combined_chunk = pd.concat(dfs_to_concat, axis=1, join="inner")
            if index == 0:
                combined_chunk.to_csv(args.output_name, mode="w", index=True)
            else:
                combined_chunk.to_csv(
                    args.output_name, mode="a", index=True, header=False
                )
