import pandas as pd
import umap

# Read beta values
beta = pd.read_csv("beta.csv", index_col=0)

# Perform UMAP
embedding = umap.UMAP().fit_transform(beta.T)

# Save UMAP embedding
umap = pd.DataFrame(data=embedding, index=beta.T.index, columns=["UMAP1", "UMAP2"])
umap.index_name = "sample"
umap.to_csv("umap.csv")