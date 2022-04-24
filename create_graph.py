import numpy as np
import scanpy as sp

cell_x_marker = np.load("CellxMarker.npy")
gene_names = np.load("GeneNames.npy")

adata = sc.AnnData(cell_x_marker)
adata.var_names = gene_names