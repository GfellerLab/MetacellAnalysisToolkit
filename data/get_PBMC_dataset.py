# Script to download pbmc3k dataset
# We need both pbmc3k and pbmc3k_processed in order to create an object containing processed cells with raw counts

import scanpy as sc 
import os

adata = sc.datasets.pbmc3k()
adata_proc = sc.datasets.pbmc3k_processed()

adata       = adata[adata_proc.obs_names].copy()
adata.obs   = adata_proc.obs.copy()
adata.uns   = adata_proc.uns.copy()
adata.obsm  = adata_proc.obsm.copy()
adata.obsp  = adata_proc.obsp.copy()

raw_ad = sc.AnnData(adata.X.copy())
raw_ad.obs_names, raw_ad.var_names = adata.obs_names, adata.var_names
adata.raw = raw_ad

directory = os.path.join("data")

if not os.path.exists(directory):
    os.makedirs(directory)
    
adata.write_h5ad(os.path.join("data", "pbmc.h5ad"))

# We keep only the file with filtered cells and raw counts
os.remove(os.path.join("data", "pbmc3k_raw.h5ad"))
os.remove(os.path.join("data", "pbmc3k_processed.h5ad"))



