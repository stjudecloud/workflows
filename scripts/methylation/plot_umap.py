import pandas as pd
import matplotlib.pyplot as plt
import argparse

def get_args():
    parser = argparse.ArgumentParser(
        description="Plot UMAP coordinates."
    )
    parser.add_argument(
        "--output-name", type=str, help="Name for output file."
    )
    parser.add_argument(
        "--umap", type=str, help="File with UMAP coordinates."
    )

    args = parser.parse_args()

    return args

if __name__ == "__main__":
    args = get_args()

    # Read UMAP embedding
    umap = pd.read_csv(args.umap, index_col=0)

    # Plot UMAP
    plt.scatter(umap["UMAP1"], umap["UMAP2"])
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.title("UMAP Embedding")
    plt.savefig(args.output_name)