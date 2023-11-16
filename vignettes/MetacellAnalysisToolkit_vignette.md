-   <a href="#load-data" id="toc-load-data">Load data</a>
-   <a href="#compute-quantitative-metrics"
    id="toc-compute-quantitative-metrics">Compute quantitative metrics</a>
    -   <a href="#purity" id="toc-purity">Purity</a>
    -   <a href="#compactness" id="toc-compactness">Compactness</a>
    -   <a href="#separation" id="toc-separation">Separation</a>
    -   <a href="#inv" id="toc-inv">INV</a>
-   <a href="#representativeness-of-metacells"
    id="toc-representativeness-of-metacells">Representativeness of
    metacells</a>

## Load data

    library(MetacellAnalysisToolkit)
    mc_data <- MetacellAnalysisToolkit::CD34_mc
    sc_data <- MetacellAnalysisToolkit::CD34_sc

## Compute quantitative metrics

### Purity

After each metacell has been annotated to the most abundant cell
category (*e.g.* cell type) composing the metacell, we can compute
metacells purity. If the annotation considered is the cell type, the
**purity** of a metacell is the proportion of the most abundant cell
type within the metacell.

    mc_data$purity <- mc_purity(membership = mc_data@misc$cell_membership$membership, annotation = sc_data$celltype)
    qc_boxplot(mc.obj = mc_data, qc.metrics = "purity")

![](/mnt/c/Aurelie/postdoc_UNIL/Metacell_review/MetacellAnalysisToolkit/vignettes/MetacellAnalysisToolkit_vignette_files/figure-markdown_strict/compute_purity-1.png)

    qc_boxplot(mc.obj = mc_data, qc.metrics = "purity", split.by = "celltype")

![](/mnt/c/Aurelie/postdoc_UNIL/Metacell_review/MetacellAnalysisToolkit/vignettes/MetacellAnalysisToolkit_vignette_files/figure-markdown_strict/compute_purity-2.png)

### Compactness

The **compactness** of a metacell is the variance of the components
within the metacell. The lower the compactness value the better.

This metric, as well as the separation metric, are computed based on a
low embedding of the single-cell data, e.g. PCA embedding which we
generate in the next chunk.

    library(Seurat)
    sc_data <- NormalizeData(sc_data, normalization.method = "LogNormalize")
    sc_data <- FindVariableFeatures(sc_data, nfeatures = 2500)
    sc_data <- ScaleData(sc_data)

    ## Centering and scaling data matrix

    sc_data <- RunPCA(sc_data, npcs = 50, verbose = F)
    sc_data <- RunUMAP(sc_data, reduction = "pca", dims = c(1:50), n.neighbors = 15, verbose = F, min.dist = 0.5)
    UMAPPlot(sc_data, group.by = "celltype", reduction = "umap")

![](/mnt/c/Aurelie/postdoc_UNIL/Metacell_review/MetacellAnalysisToolkit/vignettes/MetacellAnalysisToolkit_vignette_files/figure-markdown_strict/sc-embessing-1.png)

    membership_df <- mc_data@misc$cell_membership
    mc_data$compactness <- mc_compactness(cell.membership = membership_df, sc.obj = sc_data,
                                          sc.reduction = "pca", n.components = 30, diffusion.components = T)

    ## Computing compactness ...

    qc_boxplot(mc.obj = mc_data, qc.metrics = "compactness")

![](/mnt/c/Aurelie/postdoc_UNIL/Metacell_review/MetacellAnalysisToolkit/vignettes/MetacellAnalysisToolkit_vignette_files/figure-markdown_strict/compute_compactness-1.png)

    qc_boxplot(mc.obj = mc_data, qc.metrics = "compactness", split.by = "celltype")

![](/mnt/c/Aurelie/postdoc_UNIL/Metacell_review/MetacellAnalysisToolkit/vignettes/MetacellAnalysisToolkit_vignette_files/figure-markdown_strict/compute_compactness-2.png)

### Separation

The **separation** of a metacell is the distance to the closest metacell
\[@SEACells\]. The higher the separation value the better.

    mc_data$separation <- mc_separation(cell.membership = membership_df, sc.obj = sc_data, sc.reduction = "pca", diffusion.components = T)

    ## Computing separation ...

    qc_boxplot(mc.obj = mc_data, qc.metrics = "separation")

![](/mnt/c/Aurelie/postdoc_UNIL/Metacell_review/MetacellAnalysisToolkit/vignettes/MetacellAnalysisToolkit_vignette_files/figure-markdown_strict/compute_separation-1.png)

    qc_boxplot(mc.obj = mc_data, qc.metrics = "separation", split.by = "celltype")

![](/mnt/c/Aurelie/postdoc_UNIL/Metacell_review/MetacellAnalysisToolkit/vignettes/MetacellAnalysisToolkit_vignette_files/figure-markdown_strict/compute_separation-2.png)

### INV

The **inner normalized variance (INV)** of a metacell is the
mean-normalized variance of gene expression within the metacell. The
lower the INV value the better. Note that it is the only metric that is
latent-space independent.

    mc_data$INV <- mc_INV(cell.membership = membership_df, sc.obj = sc_data, group.label = "membership")

    ## Computing INV ...

    qc_boxplot(mc.obj = mc_data, qc.metrics = "INV")

![](/mnt/c/Aurelie/postdoc_UNIL/Metacell_review/MetacellAnalysisToolkit/vignettes/MetacellAnalysisToolkit_vignette_files/figure-markdown_strict/compute_INV-1.png)

    qc_boxplot(mc.obj = mc_data, qc.metrics = "INV", split.by = "celltype")

![](/mnt/c/Aurelie/postdoc_UNIL/Metacell_review/MetacellAnalysisToolkit/vignettes/MetacellAnalysisToolkit_vignette_files/figure-markdown_strict/compute_INV-2.png)

## Representativeness of metacells

To visualize the metacells, we can project the metacells on the
single-cell UMAP representation using the `mc_projection()` function
(adapted from the `plot.plot_2D()` from the SEACells package). A good
metacell partition should reproduce the overall structure of the
single-cell data by uniformly representing the latent space. To use this
function we need the data at the single-cell level (or at least an
low-dimensional embedding of the data) and the single-cell membership to
each the metacell.

    mc_projection(
      sc.obj = sc_data,
      mc.obj = mc_data,
      cell.membership = membership_df,
      sc.reduction = "umap",
      sc.label = "celltype", # single cells will be colored according the sc.label
      metacell.label = "celltype" # metacells cell will be colored according the metacell.label
      )

![](/mnt/c/Aurelie/postdoc_UNIL/Metacell_review/MetacellAnalysisToolkit/vignettes/MetacellAnalysisToolkit_vignette_files/figure-markdown_strict/visualize_metacells-1.png)

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

![](/mnt/c/Aurelie/postdoc_UNIL/Metacell_review/MetacellAnalysisToolkit/vignettes/MetacellAnalysisToolkit_vignette_files/figure-markdown_strict/visualize_metacells-2.png)

By default the size of the metacells dots is proportionnal to the size
of the metacells. Metacells can also be colored by a continuous variable
such as one of the QC metrics computed in the previous chunks:

    mc_projection(
      sc.obj = sc_data,
      mc.obj = mc_data,
      cell.membership = membership_df,
      sc.reduction = "umap",
      sc.label = "celltype", # single cells will be colored according the sc.label
      continuous_metric = TRUE,
      metric = "compactness"
      )

![](/mnt/c/Aurelie/postdoc_UNIL/Metacell_review/MetacellAnalysisToolkit/vignettes/MetacellAnalysisToolkit_vignette_files/figure-markdown_strict/visualize_metacells_continuous-1.png)
