import pandas as pd
import matplotlib.pyplot as plt
import sys

# Read UMAP embedding
umap = pd.read_csv(sys.argv[1], index_col=0)

# Plot UMAP
plt.scatter(umap["UMAP1"], umap["UMAP2"])
plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.title("UMAP Embedding")
plt.savefig("umap.png")  