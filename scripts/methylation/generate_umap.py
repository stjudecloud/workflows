import pandas as pd
import umap
import argparse


def get_args():
    parser = argparse.ArgumentParser(description="Generate UMAP coordinates.")
    parser.add_argument("--output-name", type=str, help="Name for output file.")
    parser.add_argument("--beta", type=str, help="File with beta values.")

    args = parser.parse_args()

    return args


if __name__ == "__main__":
    args = get_args()

    # Read beta values
    beta = pd.read_csv(args.beta, index_col=0)

    # Perform UMAP
    embedding = umap.UMAP().fit_transform(beta.T)

    # Save UMAP embedding
    umap = pd.DataFrame(data=embedding, index=beta.T.index, columns=["UMAP1", "UMAP2"])
    umap.index_name = "sample"
    umap.to_csv(args.output_name)
