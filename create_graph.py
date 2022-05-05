# reads and merges gene expression data, copy-number-variant data, and sample metadata.
# saves adata object for all samples, and adata objects for each subject (4 total). 
import numpy as np
import scanpy as sc
import pandas as pd

# Import transcriptomic data
adata = sc.read("source-data/GBM_normalized_gene_counts.txt",delimiter=" ",first_column_names=True)
adata = adata.T

# Import metadata and annotate cells with subject source information
meta = pd.read_csv("source-data/GBM_metadata.csv",delimiter=" ")
adata.obs["subject"] = pd.Categorical(meta["Sample.name"])

# Import CNV calls exported by copyKAT
adata_cnv = sc.read("source-data/CNA_matrix.csv",first_column_names=True)
adata_cnv = adata_cnv.T

# Subset transcriptomic adata object to include only subjects and genes that passed QC in copyKAT
adata = adata[adata.obs_names.isin(adata_cnv.obs_names.str[1:]),:]
adata = adata[:, adata.var_names.isin(adata_cnv.var_names)]

# Set index in adata object to match cnv data (ordered by chr position)
adata.var = adata.var.reindex(adata_cnv.var.index)

# Annoate transcriptomic adata object with CNV information
adata.obsm["CNV"] = adata_cnv.X

# Generate KNN graph using top 10 PCs
sc.tl.pca(adata, n_comps=10)
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=10, knn=True)

# Save adata object to disk
adata.write('adata.h5ad', compression="gzip")