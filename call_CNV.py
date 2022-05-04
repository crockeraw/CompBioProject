import anndata as ad
import numpy as np

adata = ad.read("adata.h5ad")

n_samples = adata_cnv.obs.shape[0]
n_copies = 2

max_copies = np.max(adata.obsm["CNV"],axis=0)
has_cnv = np.sum(adata.obsm["CNV"] > 1, axis=0)
# for each gene, n cells with the cnv
is_cnv = np.logical_and(has_cnv > .05*n_samples, max_copies > n_copies)
indices = np.where(is_cnv)[0]
adata.obsm["isCNV"] = adata.obsm["CNV"][:, indices]

adata.write('adata_CNVcalled.h5ad', compression="gzip")