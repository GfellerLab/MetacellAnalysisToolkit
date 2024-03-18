-   <a href="#integration-of-the-human-lung-cell-atlas-with-metacells"
    id="toc-integration-of-the-human-lung-cell-atlas-with-metacells">Integration
    of the human lung cell atlas with metacells</a>
    -   <a href="#setting-up-the-environment"
        id="toc-setting-up-the-environment">Setting up the environment</a>
    -   <a href="#downloading-the-data"
        id="toc-downloading-the-data">Downloading the data</a>
    -   <a href="#splitting-atlas-by-datasets"
        id="toc-splitting-atlas-by-datasets">Splitting atlas by datasets</a>
    -   <a href="#building-metacell" id="toc-building-metacell">Building
        metacell</a>
    -   <a href="#loading-metacell-objects"
        id="toc-loading-metacell-objects">Loading metacell objects</a>
    -   <a href="#merging-objects-and-basic-quality-control"
        id="toc-merging-objects-and-basic-quality-control">Merging objects and
        basic quality control</a>
    -   <a href="#unintegrated-analysis"
        id="toc-unintegrated-analysis">Unintegrated analysis</a>
    -   <a href="#seurat-integration" id="toc-seurat-integration">Seurat
        integration</a>
    -   <a href="#downstream-analysis" id="toc-downstream-analysis">Downstream
        analysis</a>
        -   <a href="#clustering" id="toc-clustering">Clustering</a>
        -   <a href="#differently-expressed-gene-deg-analysis."
            id="toc-differently-expressed-gene-deg-analysis.">Differently expressed
            gene (DEG) analysis.</a>
        -   <a href="#cell-type-abundances-analyses."
            id="toc-cell-type-abundances-analyses.">Cell type abundances
            analyses.</a>
    -   <a href="#conclusion" id="toc-conclusion">Conclusion</a>

# Integration of the human lung cell atlas with metacells

In this example we will work with the Human Cell Lung Atlas core
[HLCA](https://www.nature.com/articles/s41591-023-02327-2) gathering
around 580,000 cells from 107 individuals distributed in 166 samples.

The aim of this tutorial is to show how you can use metacells to analyze
a very large dataset using a reasonable amount of time and memory. For
this we will use here **SuperCell** via the **MATK** command line tool.

Be sure to be in the **MetacellAnalysisToolkit** environment when you
are running this Rmarkdown.

## Setting up the environment

    library(Seurat)

    ## The legacy packages maptools, rgdal, and rgeos, underpinning this package
    ## will retire shortly. Please refer to R-spatial evolution reports on
    ## https://r-spatial.org/r/2023/05/15/evolution4.html for details.
    ## This package is now running under evolution status 0

    ## Attaching SeuratObject

    library(anndata)
    library(SuperCell)
    library(ggplot2)
    wilcox.test <- "wilcox"
    if(packageVersion("Seurat") >= 5) {
      options(Seurat.object.assay.version = "v4") 
      wilcox.test <- "wilcox_limma"
      print("you are using seurat v5 with assay option v4")}

    color.celltypes  <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                          '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                          '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', "#b20000",'#E4C755', '#F7F398',
                          '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                          '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                          '#968175')

## Downloading the data

You can download the annotated data from
[cellxgene](https://cellxgene.cziscience.com/collections/6f6d381a-7701-4781-935c-db10d30de293).
Choose the `.h5ad` option after clicking on the download button for the
core atlas (3 tissues, 584,944 cells).

You can use a bash command line of this form to download the data
directly in the `./HLCA_data` directory. You will have to update the
link (obtained by clicking on download, .h5ad selection) as links are
temporary.

Please note that this may take some time as the file is quite large (~6
GB)

    #Uncomment to download the data in the ./HLCA_data/ directory after updating the link
    #mkdir -p ./HLCA_data
    #curl -o ./HLCA_data/local.h5ad "https://corpora-data-prod.s3.amazonaws.com/7bcad396-49c3-40d9-80c1-16d74e7b88bd/local.h5ad?AWSAccessKeyId=ASIATLYQ5N5XZ2V3CYXW&Signature=CI8hgXdSO2ewDXpP%2FCb7ouxW6R8%3D&x-amz-security-token=IQoJb3JpZ2luX2VjEC0aCXVzLXdlc3QtMiJGMEQCIGoOTAGVxanApGEIeRVOL%2BRK7silMZiTtgLE%2BXguyjPjAiARoOLhXmQwzwHgme2Ll0OIZK0VIrBLaH3bSbFzRzBfuSrrAwh2EAEaDDIzMTQyNjg0NjU3NSIMbCRmBRpD%2BT0U5T8%2BKsgDcLw0fAhlIgdEjdOw%2FvUOo36uXvDClcBPXmosjNUDGVIYy67gprxvikZ%2FZHqtu%2BnodejEEIIxGJw2kv0l7dcjmGgP9IFLP6WBmsGekfI7kFCkFypmZtKXqggx9stp2K3MZCrsfcEcWttsV62c690lzdiQ4UI4lUqGqXq8C7Ah1RnxfXPQJsa3YKmHs39c3mX%2BHG5Nv4rydgzhkWE7qTkGxZvqV1cLuPMz2X78zBq5GXY0HTaGvGMgAzE5OcKbqF50sxmh0pE7PGmvz1wLYN8LB6YpMbD8qCXMdP7e4uBk2yjkK23m5m%2FrMVrCWEarSh5QqrzDR347XTg%2BkVDY301ygqy3GpCTq342sTKmUZH0PRhkliGyKvakNQU4QBy6meSQORvRX1WEhn0cRYPygyD9ugK2sDqtBl0JXUlEfqSDmE%2BXGDoRFGnKiTDSvnHhVgj64h4eTUcutZFdTILwMaYGEIl1ItElCptqvYS3rmrzdvAr5nSjx%2BnK9tKt6linyh%2Bau7zc6IfQSTzZoMut%2Fw1fOuCQ%2BQmxCaEyBXzfTTrx4%2FuxyiYAkPN0vLTtSvtuklZH7O1axMTQIonnFDsnKeVnUzl3ZEgdUbxhMLL20qoGOqYBdtJOXqTiQUDX4ZH0ReubHpog%2BorDorDJ0B08Edu6k36SwuSNu6Hv8MW%2BdWFVfqs0X%2Fx74oMs8yQC8T1gSG2HrlCfLoWIBep9lA9EHq4vUBhYB4mmJ7Fsc2MdhOtof%2BzrE8b1ILxU%2Fdeliek9Aqz0uBWcfJsEu%2FlHrC1sX4P5F8nytcLxvzCTGB43mPHeqB5DZaAKC%2FY8SmSa9CJ1Njfz8n%2FIuTLv8w%3D%3D&Expires=1700662555"

First we need to specify that we will work with the
MetacellAnalysisToolkit conda environment (needed for anndata relying on
reticulate and the MATK tool).

    library(reticulate)
    conda_env <-  conda_list()[reticulate::conda_list()$name == "MetacellAnalysisToolkit","python"]

    Sys.setenv(RETICULATE_PYTHON = conda_env)

## Splitting atlas by datasets

First we will use anndata to read in backed mode (saving a lot of
memory) the whole atlas and write one h5ad file for each dataset. This
should take less than 10 minutes.

If you are limited in time feel free to process only a subset of the
dataset.

    adata <- read_h5ad("./HLCA_data/local.h5ad",backed = "r")
    adata$var_names <- adata$var$feature_name # We will use gene short name for downstream analyses
    datasets <- unique(adata$obs$dat)

    # If you are limited in time you can process on half of the datasets (uncomment th following line)
    # datasets <- datasets[1:7]


    print(dim(adata))

    lapply(datasets,FUN =  function(x) {
      dir.create(paste0("./HLCA_data/datasets/",x),recursive = T)
      adata.dataset <- AnnData(X = adata[adata$obs$dataset == x]$raw$X,
                               var = adata[adata$obs$dataset == x]$var,
                               obs = adata[adata$obs$dataset == x]$obs)
      #This will allow us to construct supervised metacell for each cell type in each sample later in the second example
      adata.dataset$obs$ann <- as.character(adata.dataset$obs$ann_level_3)
      # For cell without an annotation at the 3rd level we will use the second level of annotation
      adata.dataset$obs$ann[adata.dataset$obs$ann_level_3 == 'None'] = as.character(adata.dataset$obs$ann_level_2[adata.dataset$obs$ann_level_3 == 'None'])
      adata.dataset$obs$ann_sample <- paste0(adata.dataset$obs$ann,"_",adata.dataset$obs$sample)
      
      write_h5ad(adata.dataset,paste0("./HLCA_data/datasets/",x,"/sc_adata.h5ad"))
    }
    )

    remove(adata)
    gc()

## Building metacell

We build metacells with the MATK command line using SuperCell
(`-t SuperCell`). To facilitate downstream analysis of the donors we
build metacells for each sample in each dataset (`-a sample`). Here we
will use 2000 highly variable genes (`-f 2000`) to compute the PCA from
which we used 50 principal components (`-m 50`) to build a k = 30
(`-k 30`) nearest neighbor graph on which the metacells are identified
using a graining level of 50 (`-g 50`). We use an adata .h5ad output
format (`-s adata`) as it is faster to write and lighter to store than a
Seurat .rds object.

This step takes around 20 min with multiple cores (`-l 6`). Be aware
that parallel processing requires more memory (32 GB of memory required
for 6 cores).

If you are limited in memory you should still be able to process the
samples by reducing the number of cores (e.g. `-l 3`) or by sequentially
processing the samples (just remove the `-l`) in a slightly longer time.

    for d in ./HLCA_data/datasets/*;
    do ../cli/MATK -t SuperCell -i $d/sc_adata.h5ad -o $d -a sample -l 6 -n 50 -f 2000 -k 30 -g 50 -s adata
    done

## Loading metacell objects

We load the .h5ad objects and directly convert them in Seurat objects to
benefit from all the functions of this framework.

    metacell.files <- sapply(datasets, FUN = function(x){paste0("./HLCA_data/datasets/",x,"/mc_adata.h5ad")})

    metacell.objs <- lapply(X = metacell.files, function(X){
      adata <- read_h5ad(X)
      countMatrix <- Matrix::t(adata$X)
      colnames(countMatrix) <- adata$obs_names
      rownames(countMatrix) <- adata$var_names
      sobj <- Seurat::CreateSeuratObject(counts = countMatrix,meta.data = adata$obs)
      if(packageVersion("Seurat") >= 5){sobj[["RNA"]] <- as(object = sobj[["RNA"]], Class = "Assay")}
      sobj <- RenameCells(sobj, add.cell.id = unique(sobj$sample)) # we give unique name to metacells
      return(sobj)
    })

## Merging objects and basic quality control

Given the single-cell metadata, the MATK tool automatically assign
annotations to metacells and computes purities for all the categorical
variables present in the metadata of the input single-cell object.

Thus, we can check the purity of our metacells at different levels of
annotations, as well as their size (number of single cells they
contain).

To do so we merge the object together and use the Seurat `VlnPlot`
function.

    unintegrated.mc <- merge(metacell.objs[[1]],metacell.objs[-1])


    VlnPlot(unintegrated.mc,features = c("size","ann_level_1_purity"),group.by = 'dataset',pt.size = 0.001,ncol=2)

![](HLCA_core_atlas_files/figure-markdown_strict/unnamed-chunk-7-1.png)

    VlnPlot(unintegrated.mc,features = c("ann_level_2_purity","ann_level_3_purity"),group.by = 'dataset',pt.size = 0.001,ncol=2)

![](HLCA_core_atlas_files/figure-markdown_strict/unnamed-chunk-7-2.png)

We can also use box plots.

    ggplot(unintegrated.mc@meta.data,aes(x=dataset,y=ann_level_2_purity,fill = dataset)) + geom_boxplot() +
      scale_x_discrete(guide = guide_axis(angle = 45)) 

![](HLCA_core_atlas_files/figure-markdown_strict/unnamed-chunk-8-1.png)

    ggplot(unintegrated.mc@meta.data,aes(x=dataset,y=ann_level_3_purity,fill = dataset)) + geom_boxplot() +
      scale_x_discrete(guide = guide_axis(angle = 45)) 

![](HLCA_core_atlas_files/figure-markdown_strict/unnamed-chunk-8-2.png)

    ggplot(unintegrated.mc@meta.data,aes(x=dataset,y=ann_level_4_purity,fill = dataset)) + geom_boxplot() +
      scale_x_discrete(guide = guide_axis(angle = 45)) 

![](HLCA_core_atlas_files/figure-markdown_strict/unnamed-chunk-8-3.png)

    ggplot(unintegrated.mc@meta.data,aes(x=dataset,y=ann_finest_level_purity,fill = dataset)) + geom_boxplot() +
      scale_x_discrete(guide = guide_axis(angle = 45)) 

![](HLCA_core_atlas_files/figure-markdown_strict/unnamed-chunk-8-4.png)

Overall metacells from the different datasets present a good purity
until the third level of annotation.

## Unintegrated analysis

Let’s first do a standard dimensionality reduction without batch
correction.

    DefaultAssay(unintegrated.mc) <- "RNA"
    unintegrated.mc <- NormalizeData(unintegrated.mc)
    unintegrated.mc <- FindVariableFeatures(unintegrated.mc)
    unintegrated.mc <- ScaleData(unintegrated.mc)
    unintegrated.mc <- RunPCA(unintegrated.mc)
    unintegrated.mc <- RunUMAP(unintegrated.mc,dims = 1:30)

    umap.unintegrated.datasets <- DimPlot(unintegrated.mc,reduction = "umap",group.by = "dataset") + NoLegend() + ggtitle("unintegrated datasets")
    umap.unintegrated.types <- DimPlot(unintegrated.mc,reduction = "umap",group.by = "ann_level_2",label = T,repel = T,cols = color.celltypes)+ NoLegend() + ggtitle("unintegrated cell types")

    umap.unintegrated.datasets + umap.unintegrated.types

![](HLCA_core_atlas_files/figure-markdown_strict/unnamed-chunk-9-1.png)

    remove(unintegrated.mc) # we won't use the unintegrated object anymore
    gc()

You can see on the plots that a batch effect is clearly present at the
metacell level with metacells clustering by datasets inside the major
cell types. Let’s correct it.

## Seurat integration

Here we will use the standard Seurat\_v4 batch correction
[workflow](https://satijalab.org/seurat/archive/v4.3/integration_rpca).
As in the original study, we use the dataset rather than the donor as
the batch parameter. See method section “Data integration benchmarking”
of the [original
study](https://www.nature.com/articles/s41591-023-02327-2) for more
details.

This should take less than 5 minutes.

    n.metacells <- sapply(metacell.objs,FUN = function(x){ncol(x)})
    names(n.metacells) <- datasets
    ref.names <- sort(n.metacells,decreasing = T)[1:5]
    ref.index <- which(datasets %in% names(ref.names))

    # normalize each dataset 
    metacell.objs <- lapply(X = metacell.objs, FUN = function(x) {
      DefaultAssay(x) <- "RNA";
      x <- RenameCells(x, add.cell.id = unique(x$sample)) # we give unique name to metacells
      x <- NormalizeData(x)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
      return(x)})

    features <- SelectIntegrationFeatures(object.list = metacell.objs)

    metacell.objs <- lapply(X = metacell.objs, FUN = function(x) {
        x <- ScaleData(x, features = features, verbose = FALSE)
        x <- RunPCA(x, features = features, verbose = FALSE)
    })

    anchors <- FindIntegrationAnchors(object.list = metacell.objs, 
                                           anchor.features = features, 
                                           reduction = "rpca",
                                           reference = ref.index, # the 5 biggest datasets (in term of metacell number) are used as reference
                                           dims = 1:30)

    remove(metacell.objs) # We don't need the object list anymore
    gc()

    combined.mc <- IntegrateData(anchorset = anchors,k.weight = 40) # we have to update the k.weight parameters because the smallest dataset contain less than 100 metacells

Check the obtained object.

    combined.mc

    ## An object of class Seurat 
    ## 30024 features across 11706 samples within 2 assays 
    ## Active assay: integrated (2000 features, 2000 variable features)
    ##  1 other assay present: RNA

We can verify that the sum of metacell sizes correspond to the original
number of single-cells.

    sum(combined.mc$size)

    ## [1] 584944

Seurat returns the slot `"integrated"` that we can use for the
downstream analysis.

    DefaultAssay(combined.mc) = "integrated"
    combined.mc <- ScaleData(combined.mc, verbose = FALSE)
    combined.mc <- RunPCA(combined.mc, npcs = 30, verbose = FALSE)
    combined.mc <- RunUMAP(combined.mc, reduction = "pca", dims = 1:30)

    ## 18:00:20 UMAP embedding parameters a = 0.9922 b = 1.112

    ## 18:00:20 Read 11706 rows and found 30 numeric columns

    ## 18:00:20 Using Annoy for neighbor search, n_neighbors = 30

    ## 18:00:20 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 18:00:21 Writing NN index file to temp file /tmp/39640025/RtmpTB8dcq/file3645457b341dad
    ## 18:00:21 Searching Annoy index using 1 thread, search_k = 3000
    ## 18:00:24 Annoy recall = 100%
    ## 18:00:24 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
    ## 18:00:25 Initializing from normalized Laplacian + noise (using irlba)
    ## 18:00:26 Commencing optimization for 200 epochs, with 480248 positive edges
    ## 18:00:30 Optimization finished

    combined.mc <- RunUMAP(combined.mc, dims = 1:30,reduction =  "pca",reduction.name = "umap")

    ## 18:00:30 UMAP embedding parameters a = 0.9922 b = 1.112
    ## 18:00:30 Read 11706 rows and found 30 numeric columns
    ## 18:00:30 Using Annoy for neighbor search, n_neighbors = 30
    ## 18:00:30 Building Annoy index with metric = cosine, n_trees = 50
    ## 0%   10   20   30   40   50   60   70   80   90   100%
    ## [----|----|----|----|----|----|----|----|----|----|
    ## **************************************************|
    ## 18:00:31 Writing NN index file to temp file /tmp/39640025/RtmpTB8dcq/file36454540153459
    ## 18:00:31 Searching Annoy index using 1 thread, search_k = 3000
    ## 18:00:33 Annoy recall = 100%
    ## 18:00:34 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
    ## 18:00:34 Initializing from normalized Laplacian + noise (using irlba)
    ## 18:00:35 Commencing optimization for 200 epochs, with 480248 positive edges
    ## 18:00:39 Optimization finished

Now we can make the plots and visually compare the results with the
unintegrated analysis.

    umap.integrated.datasets <- DimPlot(combined.mc,reduction = "umap",group.by = "dataset") + NoLegend() + ggtitle("integrated datasets")
    umap.integrated.celltypes <- DimPlot(combined.mc,reduction = "umap",group.by = "ann_level_2",label = T,repel = T,cols = color.celltypes) + NoLegend() + ggtitle("integrated cell types")

    umap.integrated.datasets + umap.integrated.celltypes + umap.unintegrated.datasets + umap.unintegrated.types

![](HLCA_core_atlas_files/figure-markdown_strict/unnamed-chunk-14-1.png)

Seurat efficiently corrected the batch effect in the data while keeping
the cell type separated, but other batch correction methods such as
harmony would have also done the job.

Note that In the original study, datasets were integrated using SCANVI
semi-supervised integration using partial annotation obtained for each
dataset prior integration. If you are interested in such supervised
approach at the metacell level in R you can have a look to our second
[example](./HLCA_core_atlas_supervised.md) using the
[STACAS](https://github.com/carmonalab/STACAS) package.

We can navigate in the different annotation levels.

    library(ggplot2)

    DimPlot(combined.mc,group.by = "ann_level_1",reduction = "umap",label = T,repel = T,cols= color.celltypes) + NoLegend()

![](HLCA_core_atlas_files/figure-markdown_strict/unnamed-chunk-15-1.png)

    DimPlot(combined.mc,group.by = "ann_level_2",reduction = "umap",label = T,repel = T,cols= color.celltypes) + NoLegend()

![](HLCA_core_atlas_files/figure-markdown_strict/unnamed-chunk-15-2.png)

    DimPlot(combined.mc,group.by = "ann_level_3",reduction = "umap",label = T, repel = T,cols= color.celltypes) + NoLegend()

![](HLCA_core_atlas_files/figure-markdown_strict/unnamed-chunk-15-3.png)

## Downstream analysis

### Clustering

We cluster the metacells based on the corrected PCA space by Seurat.

    DefaultAssay(combined.mc) <- "integrated"
    combined.mc <- FindNeighbors(combined.mc,reduction = "pca",dims = 1:30)
    combined.mc <- FindClusters(combined.mc,resolution = 0.5) 
    UMAPPlot(combined.mc,label=T) + NoLegend()

![](HLCA_core_atlas_files/figure-markdown_strict/unnamed-chunk-16-1.png)

### Differently expressed gene (DEG) analysis.

We can find which of the obtained clusters correspond to the B cells
annotated in the original study (4th level of annotation)

    b.clust <- names(which.max(table(combined.mc$ann_level_4,combined.mc$integrated_snn_res.0.5)["B cells",]))
    b.clust

    ## [1] "23"

Now let’s found the markers of this cluster.

    DefaultAssay(combined.mc) <- "RNA"


    markersB<- FindMarkers(combined.mc,ident.1 = b.clust,only.pos = T,logfc.threshold = 0.25,test.use =  wilcox.test)

    ## For a more efficient implementation of the Wilcoxon Rank Sum Test,
    ## (default method for FindMarkers) please install the limma package
    ## --------------------------------------------
    ## install.packages('BiocManager')
    ## BiocManager::install('limma')
    ## --------------------------------------------
    ## After installation of limma, Seurat will automatically use the more 
    ## efficient implementation (no further action necessary).
    ## This message will be shown once per session

    head(markersB)

    ##       p_val avg_log2FC pct.1 pct.2 p_val_adj
    ## TCL1A     0  1.3939363 0.695 0.020         0
    ## FCRLA     0  1.0706003 0.962 0.022         0
    ## BLK       0  1.2209493 0.990 0.058         0
    ## FCRL5     0  0.7922142 0.886 0.027         0
    ## PNOC      0  0.5425023 0.924 0.050         0
    ## PAX5      0  0.6124840 0.895 0.023         0

This cluster clearly present a B cell signature with marker genes such
as CD19 and PAX5

    genes <-c("CD19","PAX5") # knwon mast cells markers 
    markersB[genes,]

    ##              p_val avg_log2FC pct.1 pct.2     p_val_adj
    ## CD19 6.500657e-208   1.104315 0.990 0.113 1.821744e-203
    ## PAX5  0.000000e+00   0.612484 0.895 0.023  0.000000e+00

    VlnPlot(combined.mc,genes,ncol = 1)

![](HLCA_core_atlas_files/figure-markdown_strict/unnamed-chunk-19-1.png)

By looking at the metacell annotation (assigned from the original
single-cell metadata by MATK), we can verify that we correctly retrieved
the B cell lineage cluster

    DimPlot(combined.mc[,combined.mc$integrated_snn_res.0.5 == b.clust],group.by = c("ann_level_4","integrated_snn_res.0.5"),ncol = 2)

![](HLCA_core_atlas_files/figure-markdown_strict/unnamed-chunk-20-1.png)

### Cell type abundances analyses.

We can easily make analysis of cell type abundances for different
clinical variables as we construct metacell by sample. We have to take
metacell size into account for these analyses. For instance we can
analyse the proportion of different epithelial cell types depending on
the smoking status.

    library(reshape2)
    combined.mc.epith <- combined.mc[,combined.mc$ann_level_1 == "Epithelial"] 
    #combined.metacells$major_type <- droplevels(combined.metacells$major_type)
    smpCounts <- aggregate(combined.mc.epith$size, by=list(sample = combined.mc.epith$sample,
                                                            major_type = combined.mc.epith$ann_level_3,
                                                            smoking_status = combined.mc.epith$smoking_status),
                                                            FUN=sum)


    contingencyTable <- xtabs(x ~ major_type+smoking_status,data = smpCounts)

    freqMatrix <- apply(contingencyTable,1,FUN = function(x){x/colSums(contingencyTable)})
    # res <- chisq.test(contingencyTable)
    # Roe <- res$observed/res$expected

    freqMatrix_df <- melt(freqMatrix)

    remove(combined.mc.epith)
    gc()

    ggplot(freqMatrix_df,aes(y=value,x = smoking_status,fill=major_type))  + geom_bar(stat="identity", position="fill") + scale_fill_manual(values = color.celltypes) + scale_y_continuous(labels = scales::percent) + ylab("% epithelial cells")

![](HLCA_core_atlas_files/figure-markdown_strict/unnamed-chunk-21-1.png)

Samples from smokers seem to present more AT2 cells but this quick
analysis is for illustrative purposes only. In practice it’s far more
complex to draw conclusion as we should have considered the variations
between samples/donors as well as many other technical (tissue
dissociation protocol, tissue sampling method, single-cell platform, … )
and biological (BMI, sex, Age, …) variables.

## Conclusion

Overall we made a precise simplification of the original atlas using
metacells built from each sample separately. By reducing the size of the
original atlas by a factor of 50 we could load the data, make an
integration to correct batch effect and recapitulate the main different
cell types using a reasonable amount of time and memory. In contrast,
simply loading the original single-cell data in R using Seurat is
extremely time-consuming and challenging even for the most powerful
computers.

In this first example we used a fully unsupervised workflow and did not
use any prior biological knowledge. Authors of the original study made a
remarkable work annotating the hundreds of thousands cells of the atlas.
In the second [example](./HLCA_core_atlas_supervised.md) we propose a
supervised workflow using this annotation to guide both the metacell
identification and the batch correction.

We can save the results for comparison with the second example.

    saveRDS(combined.mc,"./HLCA_data/combined.mc.unsup.rds")
