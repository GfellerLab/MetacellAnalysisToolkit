# Rscript to obtain .rds object of the pbmc dataset (to be run after get_PBMC_dataset.py)
library(reticulate)
library(Seurat)
library(anndata)
adata <- read_h5ad("data/pbmc.h5ad")
raw_counts <- Matrix::t(as(adata$raw$X, "CsparseMatrix"))
colnames(raw_counts) <- rownames(adata$obs)
rownames(raw_counts) <- rownames(adata$var)

pbmc <- CreateSeuratObject(counts = raw_counts, meta.data = adata$obs)

saveRDS(pbmc, file = paste0("data/pbmc.rds"))