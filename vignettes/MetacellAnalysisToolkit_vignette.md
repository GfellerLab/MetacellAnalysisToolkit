---
title: "MetacellAnalysisToolkit package"
author: "Aurélie Gabriel, Léonard Hérault, Mariia Bilous and David Gfeller"
date: "2023-12-20"
output:
  md_document:
    toc: true # table of content true
vignette: >
  %\VignetteIndexEntry{MetacellAnalysisToolkit package}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



## Load data


```r
library(MetacellAnalysisToolkit)
mc_data <- MetacellAnalysisToolkit::CD34_mc
sc_data <- MetacellAnalysisToolkit::CD34_sc
```

## Compute quantitative metrics

### Purity
After each metacell has been annotated to the most abundant cell category (*e.g.* cell type) composing the metacell, we can compute metacells purity. 
If the annotation considered is the cell type, the **purity** of a metacell is the proportion of the most abundant
cell type within the metacell.

```r
mc_data$purity <- mc_purity(membership = mc_data@misc$cell_membership$membership, annotation = sc_data$celltype)
qc_boxplot(mc.obj = mc_data, qc.metrics = "purity")
```

![plot of chunk compute_purity](./MetacellAnalysisToolkit_vignette_files/compute_purity-1.png)

```r
qc_boxplot(mc.obj = mc_data, qc.metrics = "purity", split.by = "celltype")
```

![plot of chunk compute_purity](./MetacellAnalysisToolkit_vignette_files/compute_purity-2.png)

### Compactness
The **compactness** of a metacell is the variance of the components within the metacell.
The lower the compactness value the better.

This metric, as well as the separation metric, are computed based on a low embedding of the single-cell data, e.g. PCA embedding 
which we generate in the next chunk.


```r
library(Seurat)
sc_data <- NormalizeData(sc_data, normalization.method = "LogNormalize")
sc_data <- FindVariableFeatures(sc_data, nfeatures = 2500)
sc_data <- ScaleData(sc_data)
```

```
## Centering and scaling data matrix
```

```r
sc_data <- RunPCA(sc_data, npcs = 50, verbose = F)
sc_data <- RunUMAP(sc_data, reduction = "pca", dims = c(1:50), n.neighbors = 15, verbose = F, min.dist = 0.5)
UMAPPlot(sc_data, group.by = "celltype", reduction = "umap")
```

![plot of chunk sc-embessing](./MetacellAnalysisToolkit_vignette_files/sc-embessing-1.png)

We compute these metrics in a diffusion map obtained from the pca.

```r
membership_df <- mc_data@misc$cell_membership
diffusion_comp <- get_diffusion_comp(sc.obj = sc_data, dims = 1:30)
```

```
## Computing diffusion maps ...
```

```
## Determing nearest neighbor graph...
```

```r
mc_data$compactness <- mc_compactness(cell.membership = membership_df,
                                      sc.obj = sc_data,
                                      sc.reduction = diffusion_comp)
qc_boxplot(mc.obj = mc_data, qc.metrics = "compactness")
```

![plot of chunk compute_compactness](./MetacellAnalysisToolkit_vignette_files/compute_compactness-1.png)

```r
qc_boxplot(mc.obj = mc_data, qc.metrics = "compactness", split.by = "celltype")
```

![plot of chunk compute_compactness](./MetacellAnalysisToolkit_vignette_files/compute_compactness-2.png)

### Separation
The **separation** of a metacell is the distance to the closest metacell [@SEACells].
The higher the separation value the better.


```r
mc_data$separation <- mc_separation(cell.membership = membership_df, 
                                    sc.obj = sc_data, 
                                    sc.reduction = diffusion_comp)
qc_boxplot(mc.obj = mc_data, qc.metrics = "separation")
```

![plot of chunk compute_separation](./MetacellAnalysisToolkit_vignette_files/compute_separation-1.png)

```r
qc_boxplot(mc.obj = mc_data, qc.metrics = "separation", split.by = "celltype")
```

![plot of chunk compute_separation](./MetacellAnalysisToolkit_vignette_files/compute_separation-2.png)


### INV
The **inner normalized variance (INV)** of a metacell is the mean-normalized variance of gene expression within the metacell.
The lower the INV value the better.
Note that it is the only metric that is latent-space independent.


```r
mc_data$INV <- mc_INV(cell.membership = membership_df, sc.obj = sc_data, group.label = "membership")
```

```
## Computing INV ...
```

```r
qc_boxplot(mc.obj = mc_data, qc.metrics = "INV")
```

![plot of chunk compute_INV](./MetacellAnalysisToolkit_vignette_files/compute_INV-1.png)

```r
qc_boxplot(mc.obj = mc_data, qc.metrics = "INV", split.by = "celltype")
```

![plot of chunk compute_INV](./MetacellAnalysisToolkit_vignette_files/compute_INV-2.png)

## Representativeness of metacells
To visualize the metacells, we can project the metacells on the single-cell UMAP representation using the `mc_projection()` function (adapted from the `plot.plot_2D()` from the SEACells package).
A good metacell partition should reproduce the overall structure of the single-cell data by uniformly representing the latent space.
To use this function we need the data at the single-cell level (or at least an low-dimensional embedding of the data) and the single-cell membership to each the metacell.


```r
mc_projection(
  sc.obj = sc_data,
  mc.obj = mc_data,
  cell.membership = membership_df,
  sc.reduction = "umap",
  sc.label = "celltype", # single cells will be colored according the sc.label
  metacell.label = "celltype" # metacells cell will be colored according the metacell.label
  )
```

![plot of chunk visualize_metacells](./MetacellAnalysisToolkit_vignette_files/visualize_metacells-1.png)

```r
# with custom colors:
colors <- c("HMP" = "#BC80BD", "DCPre" = "#66A61E", "HSC" = "#A6761D",
            "Ery" = "#E41A1C", "MEP" = "#B3B3B3", "cDC" = "#A6D854", "CLP" = "#1F78B4", "Mono" = "#E6AB02", "pDC" = "#B2DF8A" )
mc_projection(sc.obj = sc_data,
  mc.obj = mc_data,
  cell.membership = membership_df,
  sc.reduction = "umap",
  sc.label = "celltype", # single cells will be colored according the sc.label
  metacell.label = "celltype", # metacells cell will be colored according the metacell.label
  sc.color = colors,
  mc.color = colors)
```

![plot of chunk visualize_metacells](./MetacellAnalysisToolkit_vignette_files/visualize_metacells-2.png)

By default the size of the metacells dots is proportionnal to the size of the metacells.
Metacells can also be colored by a continuous variable such as one of the QC metrics computed in the previous chunks:


```r
mc_projection(
  sc.obj = sc_data,
  mc.obj = mc_data,
  cell.membership = membership_df,
  sc.reduction = "umap",
  sc.label = "celltype", # single cells will be colored according the sc.label
  continuous_metric = TRUE,
  metric = "compactness"
  )
```

![plot of chunk visualize_metacells_continuous](./MetacellAnalysisToolkit_vignette_files/visualize_metacells_continuous-1.png)
