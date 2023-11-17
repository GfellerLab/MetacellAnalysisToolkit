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
    -   <a href="#deferentially-expressed-gene-deg-analysis."
        id="toc-deferentially-expressed-gene-deg-analysis.">Deferentially
        expressed gene (DEG) analysis.</a>
    -   <a href="#cell-type-abundances-analyses."
        id="toc-cell-type-abundances-analyses.">Cell type abundances
        analyses.</a>
-   <a href="#conclusion" id="toc-conclusion">Conclusion</a>

In this example we will work with the Human Cell Lung Atlas core
[HLCA](https://www.nature.com/articles/s41591-023-02327-2) gathering
around 580,000 cells from 107 individuals distributed in 166 samples.

The aim of this tutorial is to show how you can use metacells to analyze
a very large dataset using a reasonable amount of time and memory. For
this we will use here **SuperCell** via the **MCAT** command line tool.

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
core atlas (3 tissues, 584?944 cells).

You can use a bash command line of this form to download the data
directly in the `./HLCA_data` directory. You will have to update the
link (obtained by clicking on download, .h5ad selection) as links are
temporary.

Please note that this may take some time (~45 mins) as the file is quite
large (5.6 GB)

    #Uncomment to download the data in the ./HLCA_data/ directory after updating the link
    #mkdir -p ./HLCA_data
    #curl -o ./HLCA_data/local.h5ad "https://corpora-data-prod.s3.amazonaws.com/7bcad396-49c3-40d9-80c1-16d74e7b88bd/local.h5ad?AWSAccessKeyId=ASIATLYQ5N5XZ2V3CYXW&Signature=CI8hgXdSO2ewDXpP%2FCb7ouxW6R8%3D&x-amz-security-token=IQoJb3JpZ2luX2VjEC0aCXVzLXdlc3QtMiJGMEQCIGoOTAGVxanApGEIeRVOL%2BRK7silMZiTtgLE%2BXguyjPjAiARoOLhXmQwzwHgme2Ll0OIZK0VIrBLaH3bSbFzRzBfuSrrAwh2EAEaDDIzMTQyNjg0NjU3NSIMbCRmBRpD%2BT0U5T8%2BKsgDcLw0fAhlIgdEjdOw%2FvUOo36uXvDClcBPXmosjNUDGVIYy67gprxvikZ%2FZHqtu%2BnodejEEIIxGJw2kv0l7dcjmGgP9IFLP6WBmsGekfI7kFCkFypmZtKXqggx9stp2K3MZCrsfcEcWttsV62c690lzdiQ4UI4lUqGqXq8C7Ah1RnxfXPQJsa3YKmHs39c3mX%2BHG5Nv4rydgzhkWE7qTkGxZvqV1cLuPMz2X78zBq5GXY0HTaGvGMgAzE5OcKbqF50sxmh0pE7PGmvz1wLYN8LB6YpMbD8qCXMdP7e4uBk2yjkK23m5m%2FrMVrCWEarSh5QqrzDR347XTg%2BkVDY301ygqy3GpCTq342sTKmUZH0PRhkliGyKvakNQU4QBy6meSQORvRX1WEhn0cRYPygyD9ugK2sDqtBl0JXUlEfqSDmE%2BXGDoRFGnKiTDSvnHhVgj64h4eTUcutZFdTILwMaYGEIl1ItElCptqvYS3rmrzdvAr5nSjx%2BnK9tKt6linyh%2Bau7zc6IfQSTzZoMut%2Fw1fOuCQ%2BQmxCaEyBXzfTTrx4%2FuxyiYAkPN0vLTtSvtuklZH7O1axMTQIonnFDsnKeVnUzl3ZEgdUbxhMLL20qoGOqYBdtJOXqTiQUDX4ZH0ReubHpog%2BorDorDJ0B08Edu6k36SwuSNu6Hv8MW%2BdWFVfqs0X%2Fx74oMs8yQC8T1gSG2HrlCfLoWIBep9lA9EHq4vUBhYB4mmJ7Fsc2MdhOtof%2BzrE8b1ILxU%2Fdeliek9Aqz0uBWcfJsEu%2FlHrC1sX4P5F8nytcLxvzCTGB43mPHeqB5DZaAKC%2FY8SmSa9CJ1Njfz8n%2FIuTLv8w%3D%3D&Expires=1700662555"

First we need to specify that we will work with the
MetacellAnalysisToolkit conda environment (needed for anndata relying on
reticulate and the MCAT tool).

    library(reticulate)
    conda_env <-  conda_list()[reticulate::conda_list()$name == "MetacellAnalysisToolkit","python"]

    use_condaenv(conda_env)

## Splitting atlas by datasets

First we will use anndata to read in backed mode (saving a lot of
memory) the whole atlas and write one h5ad file for each dataset. This
should take less than 10 minutes.

If you are limited in time feel free to process only a subset of the
dataset

    t0.split <- Sys.time()

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

    tf.split <- Sys.time()

    tf.split - t0.split

## Building metacell

We build metacells with the MCAT command line using SuperCell
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

    start=`date +%s`
    for d in ./HLCA_data/datasets/*;
    do ../cli/MCAT -t SuperCell -i $d/sc_adata.h5ad -o $d -a sample -l 6 -n 50 -f 2000 -k 30 -g 50 -s adata
    done
    echo "Duration: $((($(date +%s)-$start)/60)) minutes"

    ## SuperCell
    ## ./HLCA_data/datasets/Banovich_Kropski_2020/sc_adata.h5ad
    ## Identifying metacells...
    ## $ARGS
    ## character(0)
    ## 
    ## $input
    ## [1] "./HLCA_data/datasets/Banovich_Kropski_2020/sc_adata.h5ad"
    ## 
    ## $outdir
    ## [1] "./HLCA_data/datasets/Banovich_Kropski_2020"
    ## 
    ## $nPCs
    ## [1] 50
    ## 
    ## $nFeatures
    ## [1] 2000
    ## 
    ## $gamma
    ## [1] 50
    ## 
    ## $output
    ## [1] "adata"
    ## 
    ## $annotations
    ## [1] "sample"
    ## 
    ## $cores
    ## [1] 6
    ## 
    ## $k.knn
    ## [1] 30
    ## 
    ## $isNorm
    ## [1] FALSE
    ## 
    ## Loading required package: foreach
    ## Loading required package: iterators
    ## Loading required package: parallel
    ## The legacy packages maptools, rgdal, and rgeos, underpinning this package
    ## will retire shortly. Please refer to R-spatial evolution reports on
    ## https://r-spatial.org/r/2023/05/15/evolution4.html for details.
    ## This package is now running under evolution status 0 
    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
    ##             used   (Mb) gc trigger   (Mb)  max used   (Mb)
    ## Ncells   3161862  168.9    6273702  335.1   6273702  335.1
    ## Vcells 281732323 2149.5  944483300 7205.9 755448546 5763.7
    ## Normalize data...NULL
    ## Identify Metacells...            used   (Mb) gc trigger    (Mb)   max used   (Mb)
    ## Ncells   3190314  170.4    6273702   335.1    6273702  335.1
    ## Vcells 575057161 4387.4 1360231952 10377.8 1141961298 8712.5
    ## Assign metadata to metacells and compute purities...Done.
    ## SuperCell
    ## ./HLCA_data/datasets/Barbry_Leroy_2020/sc_adata.h5ad
    ## Identifying metacells...
    ## $ARGS
    ## character(0)
    ## 
    ## $input
    ## [1] "./HLCA_data/datasets/Barbry_Leroy_2020/sc_adata.h5ad"
    ## 
    ## $outdir
    ## [1] "./HLCA_data/datasets/Barbry_Leroy_2020"
    ## 
    ## $nPCs
    ## [1] 50
    ## 
    ## $nFeatures
    ## [1] 2000
    ## 
    ## $gamma
    ## [1] 50
    ## 
    ## $output
    ## [1] "adata"
    ## 
    ## $annotations
    ## [1] "sample"
    ## 
    ## $cores
    ## [1] 6
    ## 
    ## $k.knn
    ## [1] 30
    ## 
    ## $isNorm
    ## [1] FALSE
    ## 
    ## Loading required package: foreach
    ## Loading required package: iterators
    ## Loading required package: parallel
    ## The legacy packages maptools, rgdal, and rgeos, underpinning this package
    ## will retire shortly. Please refer to R-spatial evolution reports on
    ## https://r-spatial.org/r/2023/05/15/evolution4.html for details.
    ## This package is now running under evolution status 0 
    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
    ##             used   (Mb) gc trigger   (Mb)  max used   (Mb)
    ## Ncells   3113366  166.3    6204411  331.4   6204411  331.4
    ## Vcells 182491490 1392.3  609000860 4646.4 483622058 3689.8
    ## Normalize data...NULL
    ## Identify Metacells...            used   (Mb) gc trigger   (Mb)  max used   (Mb)
    ## Ncells   3140852  167.8    6204411  331.4   6204411  331.4
    ## Vcells 370466888 2826.5  877137238 6692.1 760126689 5799.4
    ## Assign metadata to metacells and compute purities...Done.
    ## SuperCell
    ## ./HLCA_data/datasets/Jain_Misharin_2021_10Xv1/sc_adata.h5ad
    ## Identifying metacells...
    ## $ARGS
    ## character(0)
    ## 
    ## $input
    ## [1] "./HLCA_data/datasets/Jain_Misharin_2021_10Xv1/sc_adata.h5ad"
    ## 
    ## $outdir
    ## [1] "./HLCA_data/datasets/Jain_Misharin_2021_10Xv1"
    ## 
    ## $nPCs
    ## [1] 50
    ## 
    ## $nFeatures
    ## [1] 2000
    ## 
    ## $gamma
    ## [1] 50
    ## 
    ## $output
    ## [1] "adata"
    ## 
    ## $annotations
    ## [1] "sample"
    ## 
    ## $cores
    ## [1] 6
    ## 
    ## $k.knn
    ## [1] 30
    ## 
    ## $isNorm
    ## [1] FALSE
    ## 
    ## Loading required package: foreach
    ## Loading required package: iterators
    ## Loading required package: parallel
    ## The legacy packages maptools, rgdal, and rgeos, underpinning this package
    ## will retire shortly. Please refer to R-spatial evolution reports on
    ## https://r-spatial.org/r/2023/05/15/evolution4.html for details.
    ## This package is now running under evolution status 0 
    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
    ##            used  (Mb) gc trigger  (Mb) max used  (Mb)
    ## Ncells  3051224 163.0    5114078 273.2  5114078 273.2
    ## Vcells 36979417 282.2   98207691 749.3 91437915 697.7
    ## Normalize data...NULL
    ## Identify Metacells...           used  (Mb) gc trigger   (Mb)  max used   (Mb)
    ## Ncells  3074906 164.3    5114078  273.2   5114078  273.2
    ## Vcells 70570589 538.5  145695173 1111.6 142055304 1083.8
    ## Assign metadata to metacells and compute purities...Done.
    ## SuperCell
    ## ./HLCA_data/datasets/Jain_Misharin_2021_10Xv2/sc_adata.h5ad
    ## Identifying metacells...
    ## $ARGS
    ## character(0)
    ## 
    ## $input
    ## [1] "./HLCA_data/datasets/Jain_Misharin_2021_10Xv2/sc_adata.h5ad"
    ## 
    ## $outdir
    ## [1] "./HLCA_data/datasets/Jain_Misharin_2021_10Xv2"
    ## 
    ## $nPCs
    ## [1] 50
    ## 
    ## $nFeatures
    ## [1] 2000
    ## 
    ## $gamma
    ## [1] 50
    ## 
    ## $output
    ## [1] "adata"
    ## 
    ## $annotations
    ## [1] "sample"
    ## 
    ## $cores
    ## [1] 6
    ## 
    ## $k.knn
    ## [1] 30
    ## 
    ## $isNorm
    ## [1] FALSE
    ## 
    ## Loading required package: foreach
    ## Loading required package: iterators
    ## Loading required package: parallel
    ## The legacy packages maptools, rgdal, and rgeos, underpinning this package
    ## will retire shortly. Please refer to R-spatial evolution reports on
    ## https://r-spatial.org/r/2023/05/15/evolution4.html for details.
    ## This package is now running under evolution status 0 
    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
    ##            used  (Mb) gc trigger   (Mb)  max used   (Mb)
    ## Ncells  3072019 164.1    5989752  319.9   5989752  319.9
    ## Vcells 70366982 536.9  191239283 1459.1 182966751 1396.0
    ## Normalize data...NULL
    ## Identify Metacells...            used   (Mb) gc trigger   (Mb)  max used   (Mb)
    ## Ncells   3095827  165.4    5989752  319.9   5989752  319.9
    ## Vcells 140370717 1071.0  286831697 2188.4 286051531 2182.4
    ## Assign metadata to metacells and compute purities...Done.
    ## SuperCell
    ## ./HLCA_data/datasets/Krasnow_2020/sc_adata.h5ad
    ## Identifying metacells...
    ## $ARGS
    ## character(0)
    ## 
    ## $input
    ## [1] "./HLCA_data/datasets/Krasnow_2020/sc_adata.h5ad"
    ## 
    ## $outdir
    ## [1] "./HLCA_data/datasets/Krasnow_2020"
    ## 
    ## $nPCs
    ## [1] 50
    ## 
    ## $nFeatures
    ## [1] 2000
    ## 
    ## $gamma
    ## [1] 50
    ## 
    ## $output
    ## [1] "adata"
    ## 
    ## $annotations
    ## [1] "sample"
    ## 
    ## $cores
    ## [1] 6
    ## 
    ## $k.knn
    ## [1] 30
    ## 
    ## $isNorm
    ## [1] FALSE
    ## 
    ## Loading required package: foreach
    ## Loading required package: iterators
    ## Loading required package: parallel
    ## The legacy packages maptools, rgdal, and rgeos, underpinning this package
    ## will retire shortly. Please refer to R-spatial evolution reports on
    ## https://r-spatial.org/r/2023/05/15/evolution4.html for details.
    ## This package is now running under evolution status 0 
    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
    ##             used   (Mb) gc trigger   (Mb)  max used   (Mb)
    ## Ncells   3099619  165.6    6184490  330.3   6184490  330.3
    ## Vcells 193318403 1475.0  649750587 4957.3 511849336 3905.2
    ## Normalize data...NULL
    ## Identify Metacells...            used   (Mb) gc trigger   (Mb)  max used   (Mb)
    ## Ncells   3123478  166.9    6184490  330.3   6184490  330.3
    ## Vcells 389842971 2974.3  942659982 7192.0 785461369 5992.6
    ## Assign metadata to metacells and compute purities...Done.
    ## SuperCell
    ## ./HLCA_data/datasets/Lafyatis_Rojas_2019_10Xv1/sc_adata.h5ad
    ## Identifying metacells...
    ## $ARGS
    ## character(0)
    ## 
    ## $input
    ## [1] "./HLCA_data/datasets/Lafyatis_Rojas_2019_10Xv1/sc_adata.h5ad"
    ## 
    ## $outdir
    ## [1] "./HLCA_data/datasets/Lafyatis_Rojas_2019_10Xv1"
    ## 
    ## $nPCs
    ## [1] 50
    ## 
    ## $nFeatures
    ## [1] 2000
    ## 
    ## $gamma
    ## [1] 50
    ## 
    ## $output
    ## [1] "adata"
    ## 
    ## $annotations
    ## [1] "sample"
    ## 
    ## $cores
    ## [1] 6
    ## 
    ## $k.knn
    ## [1] 30
    ## 
    ## $isNorm
    ## [1] FALSE
    ## 
    ## Loading required package: foreach
    ## Loading required package: iterators
    ## Loading required package: parallel
    ## The legacy packages maptools, rgdal, and rgeos, underpinning this package
    ## will retire shortly. Please refer to R-spatial evolution reports on
    ## https://r-spatial.org/r/2023/05/15/evolution4.html for details.
    ## This package is now running under evolution status 0 
    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
    ##           used  (Mb) gc trigger  (Mb) max used  (Mb)
    ## Ncells 3041821 162.5    5086878 271.7  5086878 271.7
    ## Vcells 9221887  70.4   19726290 150.5 16276218 124.2
    ## Normalize data...NULL
    ## Identify Metacells...           used  (Mb) gc trigger  (Mb) max used  (Mb)
    ## Ncells  3065374 163.8    5086878 271.7  5086878 271.7
    ## Vcells 13767335 105.1   24193881 184.6 24098230 183.9
    ## Assign metadata to metacells and compute purities...Done.
    ## SuperCell
    ## ./HLCA_data/datasets/Lafyatis_Rojas_2019_10Xv2/sc_adata.h5ad
    ## Identifying metacells...
    ## $ARGS
    ## character(0)
    ## 
    ## $input
    ## [1] "./HLCA_data/datasets/Lafyatis_Rojas_2019_10Xv2/sc_adata.h5ad"
    ## 
    ## $outdir
    ## [1] "./HLCA_data/datasets/Lafyatis_Rojas_2019_10Xv2"
    ## 
    ## $nPCs
    ## [1] 50
    ## 
    ## $nFeatures
    ## [1] 2000
    ## 
    ## $gamma
    ## [1] 50
    ## 
    ## $output
    ## [1] "adata"
    ## 
    ## $annotations
    ## [1] "sample"
    ## 
    ## $cores
    ## [1] 6
    ## 
    ## $k.knn
    ## [1] 30
    ## 
    ## $isNorm
    ## [1] FALSE
    ## 
    ## Loading required package: foreach
    ## Loading required package: iterators
    ## Loading required package: parallel
    ## The legacy packages maptools, rgdal, and rgeos, underpinning this package
    ## will retire shortly. Please refer to R-spatial evolution reports on
    ## https://r-spatial.org/r/2023/05/15/evolution4.html for details.
    ## This package is now running under evolution status 0 
    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
    ##            used  (Mb) gc trigger   (Mb)  max used   (Mb)
    ## Ncells  3060338 163.5    5139157  274.5   5139157  274.5
    ## Vcells 58271718 444.6  158274579 1207.6 149341027 1139.4
    ## Normalize data...NULL
    ## Identify Metacells...            used  (Mb) gc trigger   (Mb)  max used   (Mb)
    ## Ncells   3084146 164.8    5139157  274.5   5139157  274.5
    ## Vcells 114332280 872.3  283801305 2165.3 236316892 1803.0
    ## Assign metadata to metacells and compute purities...Done.
    ## SuperCell
    ## ./HLCA_data/datasets/Meyer_2019/sc_adata.h5ad
    ## Identifying metacells...
    ## $ARGS
    ## character(0)
    ## 
    ## $input
    ## [1] "./HLCA_data/datasets/Meyer_2019/sc_adata.h5ad"
    ## 
    ## $outdir
    ## [1] "./HLCA_data/datasets/Meyer_2019"
    ## 
    ## $nPCs
    ## [1] 50
    ## 
    ## $nFeatures
    ## [1] 2000
    ## 
    ## $gamma
    ## [1] 50
    ## 
    ## $output
    ## [1] "adata"
    ## 
    ## $annotations
    ## [1] "sample"
    ## 
    ## $cores
    ## [1] 6
    ## 
    ## $k.knn
    ## [1] 30
    ## 
    ## $isNorm
    ## [1] FALSE
    ## 
    ## Loading required package: foreach
    ## Loading required package: iterators
    ## Loading required package: parallel
    ## The legacy packages maptools, rgdal, and rgeos, underpinning this package
    ## will retire shortly. Please refer to R-spatial evolution reports on
    ## https://r-spatial.org/r/2023/05/15/evolution4.html for details.
    ## This package is now running under evolution status 0 
    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
    ##            used  (Mb) gc trigger   (Mb)  max used   (Mb)
    ## Ncells  3074779 164.3    5993138  320.1   5993138  320.1
    ## Vcells 85572482 652.9  234274769 1787.4 223502648 1705.2
    ## Normalize data...NULL
    ## Identify Metacells...            used   (Mb) gc trigger   (Mb)  max used   (Mb)
    ## Ncells   3099490  165.6    5993138  320.1   5993138  320.1
    ## Vcells 171036610 1305.0  351766244 2683.8 350302459 2672.6
    ## Assign metadata to metacells and compute purities...Done.
    ## SuperCell
    ## ./HLCA_data/datasets/Misharin_2021/sc_adata.h5ad
    ## Identifying metacells...
    ## $ARGS
    ## character(0)
    ## 
    ## $input
    ## [1] "./HLCA_data/datasets/Misharin_2021/sc_adata.h5ad"
    ## 
    ## $outdir
    ## [1] "./HLCA_data/datasets/Misharin_2021"
    ## 
    ## $nPCs
    ## [1] 50
    ## 
    ## $nFeatures
    ## [1] 2000
    ## 
    ## $gamma
    ## [1] 50
    ## 
    ## $output
    ## [1] "adata"
    ## 
    ## $annotations
    ## [1] "sample"
    ## 
    ## $cores
    ## [1] 6
    ## 
    ## $k.knn
    ## [1] 30
    ## 
    ## $isNorm
    ## [1] FALSE
    ## 
    ## Loading required package: foreach
    ## Loading required package: iterators
    ## Loading required package: parallel
    ## The legacy packages maptools, rgdal, and rgeos, underpinning this package
    ## will retire shortly. Please refer to R-spatial evolution reports on
    ## https://r-spatial.org/r/2023/05/15/evolution4.html for details.
    ## This package is now running under evolution status 0 
    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
    ##             used   (Mb) gc trigger   (Mb)  max used   (Mb)
    ## Ncells   3104104  165.8    6190188  330.6   6190188  330.6
    ## Vcells 253549539 1934.5  685799463 5232.3 674669700 5147.4
    ## Normalize data...NULL
    ## Identify Metacells...            used   (Mb) gc trigger   (Mb)   max used   (Mb)
    ## Ncells   3128557  167.1    6190188  330.6    6190188  330.6
    ## Vcells 511193282 3900.1 1073724162 8191.9 1071851411 8177.6
    ## Assign metadata to metacells and compute purities...Done.
    ## SuperCell
    ## ./HLCA_data/datasets/Misharin_Budinger_2018/sc_adata.h5ad
    ## Identifying metacells...
    ## $ARGS
    ## character(0)
    ## 
    ## $input
    ## [1] "./HLCA_data/datasets/Misharin_Budinger_2018/sc_adata.h5ad"
    ## 
    ## $outdir
    ## [1] "./HLCA_data/datasets/Misharin_Budinger_2018"
    ## 
    ## $nPCs
    ## [1] 50
    ## 
    ## $nFeatures
    ## [1] 2000
    ## 
    ## $gamma
    ## [1] 50
    ## 
    ## $output
    ## [1] "adata"
    ## 
    ## $annotations
    ## [1] "sample"
    ## 
    ## $cores
    ## [1] 6
    ## 
    ## $k.knn
    ## [1] 30
    ## 
    ## $isNorm
    ## [1] FALSE
    ## 
    ## Loading required package: foreach
    ## Loading required package: iterators
    ## Loading required package: parallel
    ## The legacy packages maptools, rgdal, and rgeos, underpinning this package
    ## will retire shortly. Please refer to R-spatial evolution reports on
    ## https://r-spatial.org/r/2023/05/15/evolution4.html for details.
    ## This package is now running under evolution status 0 
    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
    ##             used  (Mb) gc trigger   (Mb)  max used   (Mb)
    ## Ncells   3080327 164.6    6143747  328.2   6143747  328.2
    ## Vcells 129386820 987.2  360068441 2747.2 341284540 2603.8
    ## Normalize data...NULL
    ## Identify Metacells...            used   (Mb) gc trigger   (Mb)  max used   (Mb)
    ## Ncells   3104510  165.8    6143747  328.2   6143747  328.2
    ## Vcells 259763687 1981.9  540659870 4125.0 538085462 4105.3
    ## Assign metadata to metacells and compute purities...Done.
    ## SuperCell
    ## ./HLCA_data/datasets/Nawijn_2021/sc_adata.h5ad
    ## Identifying metacells...
    ## $ARGS
    ## character(0)
    ## 
    ## $input
    ## [1] "./HLCA_data/datasets/Nawijn_2021/sc_adata.h5ad"
    ## 
    ## $outdir
    ## [1] "./HLCA_data/datasets/Nawijn_2021"
    ## 
    ## $nPCs
    ## [1] 50
    ## 
    ## $nFeatures
    ## [1] 2000
    ## 
    ## $gamma
    ## [1] 50
    ## 
    ## $output
    ## [1] "adata"
    ## 
    ## $annotations
    ## [1] "sample"
    ## 
    ## $cores
    ## [1] 6
    ## 
    ## $k.knn
    ## [1] 30
    ## 
    ## $isNorm
    ## [1] FALSE
    ## 
    ## Loading required package: foreach
    ## Loading required package: iterators
    ## Loading required package: parallel
    ## The legacy packages maptools, rgdal, and rgeos, underpinning this package
    ## will retire shortly. Please refer to R-spatial evolution reports on
    ## https://r-spatial.org/r/2023/05/15/evolution4.html for details.
    ## This package is now running under evolution status 0 
    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
    ##             used   (Mb) gc trigger   (Mb)  max used   (Mb)
    ## Ncells   3109656  166.1    6198112  331.1   6198112  331.1
    ## Vcells 205652472 1569.1  690116576 5265.2 544745103 4156.1
    ## Normalize data...NULL
    ## Identify Metacells...            used   (Mb) gc trigger   (Mb)  max used   (Mb)
    ## Ncells   3134757  167.5    6198112  331.1   6198112  331.1
    ## Vcells 415952854 3173.5  993943869 7583.2 828080878 6317.8
    ## Assign metadata to metacells and compute purities...Done.
    ## SuperCell
    ## ./HLCA_data/datasets/Seibold_2020_10Xv2/sc_adata.h5ad
    ## Identifying metacells...
    ## $ARGS
    ## character(0)
    ## 
    ## $input
    ## [1] "./HLCA_data/datasets/Seibold_2020_10Xv2/sc_adata.h5ad"
    ## 
    ## $outdir
    ## [1] "./HLCA_data/datasets/Seibold_2020_10Xv2"
    ## 
    ## $nPCs
    ## [1] 50
    ## 
    ## $nFeatures
    ## [1] 2000
    ## 
    ## $gamma
    ## [1] 50
    ## 
    ## $output
    ## [1] "adata"
    ## 
    ## $annotations
    ## [1] "sample"
    ## 
    ## $cores
    ## [1] 6
    ## 
    ## $k.knn
    ## [1] 30
    ## 
    ## $isNorm
    ## [1] FALSE
    ## 
    ## Loading required package: foreach
    ## Loading required package: iterators
    ## Loading required package: parallel
    ## The legacy packages maptools, rgdal, and rgeos, underpinning this package
    ## will retire shortly. Please refer to R-spatial evolution reports on
    ## https://r-spatial.org/r/2023/05/15/evolution4.html for details.
    ## This package is now running under evolution status 0 
    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
    ##            used  (Mb) gc trigger   (Mb)  max used   (Mb)
    ## Ncells  3051353 163.0    5113128  273.1   5113128  273.1
    ## Vcells 70000490 534.1  193347990 1475.2 190393252 1452.6
    ## Normalize data...NULL
    ## Identify Metacells...            used   (Mb) gc trigger   (Mb)  max used   (Mb)
    ## Ncells   3074909  164.3    5113128  273.1   5113128  273.1
    ## Vcells 136566510 1042.0  288384514 2200.2 285489450 2178.2
    ## Assign metadata to metacells and compute purities...Done.
    ## SuperCell
    ## ./HLCA_data/datasets/Seibold_2020_10Xv3/sc_adata.h5ad
    ## Identifying metacells...
    ## $ARGS
    ## character(0)
    ## 
    ## $input
    ## [1] "./HLCA_data/datasets/Seibold_2020_10Xv3/sc_adata.h5ad"
    ## 
    ## $outdir
    ## [1] "./HLCA_data/datasets/Seibold_2020_10Xv3"
    ## 
    ## $nPCs
    ## [1] 50
    ## 
    ## $nFeatures
    ## [1] 2000
    ## 
    ## $gamma
    ## [1] 50
    ## 
    ## $output
    ## [1] "adata"
    ## 
    ## $annotations
    ## [1] "sample"
    ## 
    ## $cores
    ## [1] 6
    ## 
    ## $k.knn
    ## [1] 30
    ## 
    ## $isNorm
    ## [1] FALSE
    ## 
    ## Loading required package: foreach
    ## Loading required package: iterators
    ## Loading required package: parallel
    ## The legacy packages maptools, rgdal, and rgeos, underpinning this package
    ## will retire shortly. Please refer to R-spatial evolution reports on
    ## https://r-spatial.org/r/2023/05/15/evolution4.html for details.
    ## This package is now running under evolution status 0 
    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
    ##             used   (Mb) gc trigger   (Mb)  max used   (Mb)
    ## Ncells   3060346  163.5    5139852  274.5   5139852  274.5
    ## Vcells 188474317 1438.0  533042124 4066.8 496618122 3788.9
    ## Normalize data...NULL
    ## Identify Metacells...            used   (Mb) gc trigger   (Mb)  max used   (Mb)
    ## Ncells   3084541  164.8    5139852  274.5   5139852  274.5
    ## Vcells 374716173 2858.9  798693980 6093.6 765242410 5838.4
    ## Assign metadata to metacells and compute purities...Done.
    ## SuperCell
    ## ./HLCA_data/datasets/Teichmann_Meyer_2019/sc_adata.h5ad
    ## Identifying metacells...
    ## $ARGS
    ## character(0)
    ## 
    ## $input
    ## [1] "./HLCA_data/datasets/Teichmann_Meyer_2019/sc_adata.h5ad"
    ## 
    ## $outdir
    ## [1] "./HLCA_data/datasets/Teichmann_Meyer_2019"
    ## 
    ## $nPCs
    ## [1] 50
    ## 
    ## $nFeatures
    ## [1] 2000
    ## 
    ## $gamma
    ## [1] 50
    ## 
    ## $output
    ## [1] "adata"
    ## 
    ## $annotations
    ## [1] "sample"
    ## 
    ## $cores
    ## [1] 6
    ## 
    ## $k.knn
    ## [1] 30
    ## 
    ## $isNorm
    ## [1] FALSE
    ## 
    ## Loading required package: foreach
    ## Loading required package: iterators
    ## Loading required package: parallel
    ## The legacy packages maptools, rgdal, and rgeos, underpinning this package
    ## will retire shortly. Please refer to R-spatial evolution reports on
    ## https://r-spatial.org/r/2023/05/15/evolution4.html for details.
    ## This package is now running under evolution status 0 
    ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
    ##            used  (Mb) gc trigger   (Mb)  max used  (Mb)
    ## Ncells  3050792 163.0    5113418  273.1   5113418 273.1
    ## Vcells 48535686 370.3  131515210 1003.4 121596124 927.8
    ## Normalize data...NULL
    ## Identify Metacells...           used  (Mb) gc trigger   (Mb)  max used   (Mb)
    ## Ncells  3075427 164.3    5113418  273.1   5113418  273.1
    ## Vcells 93617248 714.3  195643553 1492.7 193590634 1477.0
    ## Assign metadata to metacells and compute purities...Done.
    ## Duration: 17 minutes

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
      sobj <- RenameCells(sobj, add.cell.id = unique(sobj$sample)) # we give unique name to metacells
      return(sobj)
    })

## Merging objects and basic quality control

Given the single-cell metadata, the MCAT tool automatically assign
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

Here we will use the standard Seurat\_v4 batch correction [workflow]()
As in the original study, we use the dataset rather than the donor as
the batch parameter. See method section “Data integration benchmarking”
of the [original
study](https://www.nature.com/articles/s41591-023-02327-2) for more
details.

This should take less than 5 minutes.

    # Install package if needed
    # if (!requireNamespace("STACAS")) remotes::install_github("carmonalab/STACAS")
    # 
    # library(STACAS)

    t0_integration <- Sys.time()

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
                                           reference = c(1,2,5,9,11), # the 5 biggest datasets (in term of metacell number) are used as reference
                                           dims = 1:30)

    remove(metacell.objs) # We don't need the object list anymore
    gc()

    combined.mc <- IntegrateData(anchorset = anchors,k.weight = 50) # we have to update the k.weight parameters because the smallest dataset contain less than 100 metacells



    tf_integration <- Sys.time()

    tf_integration - t0_integration

Check the obtained object.

    combined.mc

    ## An object of class Seurat 
    ## 30024 features across 11706 samples within 2 assays 
    ## Active assay: integrated (2000 features, 2000 variable features)
    ##  1 other assay present: RNA

We can verify that the sum of metacell sizes correspond to the original
number of single-cells

    sum(combined.mc$size)

    ## [1] 584944

Seurat returns the slot `"integrated"` that we can use for the
downstream analysis.

    DefaultAssay(combined.mc) = "integrated"
    combined.mc <- ScaleData(combined.mc, verbose = FALSE)
    combined.mc <- RunPCA(combined.mc, npcs = 30, verbose = FALSE)
    combined.mc <- RunUMAP(combined.mc, reduction = "pca", dims = 1:30)

    ## 17:18:05 UMAP embedding parameters a = 0.9922 b = 1.112

    ## 17:18:05 Read 11706 rows and found 30 numeric columns

    ## 17:18:05 Using Annoy for neighbor search, n_neighbors = 30

    ## 17:18:05 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 17:18:07 Writing NN index file to temp file /tmp/35267733/RtmpcpF1U0/file2301e270dfeb8
    ## 17:18:07 Searching Annoy index using 1 thread, search_k = 3000
    ## 17:18:10 Annoy recall = 100%
    ## 17:18:10 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
    ## 17:18:11 Initializing from normalized Laplacian + noise (using irlba)
    ## 17:18:12 Commencing optimization for 200 epochs, with 484456 positive edges
    ## 17:18:18 Optimization finished

    combined.mc <- RunUMAP(combined.mc, dims = 1:30,reduction =  "pca",reduction.name = "umap")

    ## 17:18:18 UMAP embedding parameters a = 0.9922 b = 1.112
    ## 17:18:18 Read 11706 rows and found 30 numeric columns
    ## 17:18:18 Using Annoy for neighbor search, n_neighbors = 30
    ## 17:18:18 Building Annoy index with metric = cosine, n_trees = 50
    ## 0%   10   20   30   40   50   60   70   80   90   100%
    ## [----|----|----|----|----|----|----|----|----|----|
    ## **************************************************|
    ## 17:18:19 Writing NN index file to temp file /tmp/35267733/RtmpcpF1U0/file2301e22d1f31c0
    ## 17:18:19 Searching Annoy index using 1 thread, search_k = 3000
    ## 17:18:22 Annoy recall = 100%
    ## 17:18:22 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
    ## 17:18:23 Initializing from normalized Laplacian + noise (using irlba)
    ## 17:18:25 Commencing optimization for 200 epochs, with 484456 positive edges
    ## 17:18:30 Optimization finished

Now we can make the plots and visually compare the results with the
unintegrated analysis.

    umap.stacas.datasets <- DimPlot(combined.mc,reduction = "umap",group.by = "dataset") + NoLegend() + ggtitle("integrated datasets")
    umap.stacas.celltypes <- DimPlot(combined.mc,reduction = "umap",group.by = "ann_level_2",label = T,repel = T,cols = color.celltypes) + NoLegend() + ggtitle("integrated cell types")

    umap.stacas.datasets + umap.stacas.celltypes + umap.unintegrated.datasets + umap.unintegrated.types

![](HLCA_core_atlas_files/figure-markdown_strict/unnamed-chunk-14-1.png)

STACAS efficiently corrected the batch effect in the data while keeping
the cell type separated, but other batch correction methods such as
harmony would have also done the job.

Note that In the original study, datasets were integrated using SCANVI
semi-supervised integration using partial annotation obtained for each
dataset prior integration. If you are interested in such supervised
approach at the metacell level in R you can have a look to our second
[example](./HLCA_core_atlas_supervised.Rmd) using the
[STACAS](https://github.com/carmonalab/STACAS) package.

We can navigate in the different annotation levels.

    library(ggplot2)

    DimPlot(combined.mc,group.by = "ann_level_1",reduction = "umap",cols= color.celltypes)

![](HLCA_core_atlas_files/figure-markdown_strict/unnamed-chunk-15-1.png)

    DimPlot(combined.mc,group.by = "ann_level_2",reduction = "umap",label = T,repel = T,cols= color.celltypes)

![](HLCA_core_atlas_files/figure-markdown_strict/unnamed-chunk-15-2.png)

    DimPlot(combined.mc,group.by = "ann_level_3",reduction = "umap",label = T, repel = T,cols= color.celltypes) + NoLegend()

![](HLCA_core_atlas_files/figure-markdown_strict/unnamed-chunk-15-3.png)

## Downstream analysis

### Clustering

We cluster the metacells based on the corrected PCA space by STACAS

    DefaultAssay(combined.mc) <- "integrated"
    combined.mc <- FindNeighbors(combined.mc,reduction = "pca",dims = 1:30)
    combined.mc <- FindClusters(combined.mc,resolution = 0.5) 
    UMAPPlot(combined.mc,label=T) + NoLegend()

![](HLCA_core_atlas_files/figure-markdown_strict/unnamed-chunk-16-1.png)

### Deferentially expressed gene (DEG) analysis.

Now let’s found the markers of the cluster 19 we’ve just identified.

    DefaultAssay(combined.mc) <- "RNA"
    markers18 <- FindMarkers(combined.mc,ident.1 = 18,only.pos = T)

    ## For a more efficient implementation of the Wilcoxon Rank Sum Test,
    ## (default method for FindMarkers) please install the limma package
    ## --------------------------------------------
    ## install.packages('BiocManager')
    ## BiocManager::install('limma')
    ## --------------------------------------------
    ## After installation of limma, Seurat will automatically use the more 
    ## efficient implementation (no further action necessary).
    ## This message will be shown once per session

    head(markers18)

    ##          p_val avg_log2FC pct.1 pct.2 p_val_adj
    ## TNFRSF17     0  0.9274364 0.762 0.029         0
    ## TCL1A        0  1.0045973 0.481 0.019         0
    ## CD79A        0  2.9087361 1.000 0.133         0
    ## FCRLA        0  0.7458042 0.801 0.018         0
    ## BLK          0  0.8766954 0.845 0.054         0
    ## FCRL5        0  0.8790523 0.934 0.020         0

This cluster clearly present a B cell signature with marker genes such
as CD19 and PAX5

    genes <-c("CD19","PAX5") # knwon mast cells markers 
    markers18[genes,]

    ##              p_val avg_log2FC pct.1 pct.2     p_val_adj
    ## CD19 3.323011e-222  0.7754585 0.812 0.111 9.312405e-218
    ## PAX5  0.000000e+00  0.3968396 0.630 0.022  0.000000e+00

    VlnPlot(combined.mc,genes,ncol = 1)

![](HLCA_core_atlas_files/figure-markdown_strict/unnamed-chunk-18-1.png)

By looking at the metacell annotation (assigned from the original
single-cell metadata by MCAT), we can verify that we correctly retrieved
the B cell lineage cluster

    DimPlot(combined.mc[,combined.mc$integrated_snn_res.0.5 == 18],group.by = c("ann_level_3","integrated_snn_res.0.5"),ncol = 2)

![](HLCA_core_atlas_files/figure-markdown_strict/unnamed-chunk-19-1.png)

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

    remove(combined.mc.epith)
    gc()

    ggplot(smpCounts,aes(x = smoking_status,fill=major_type)) + geom_bar(position = "fill") + scale_fill_manual(values = color.celltypes) + xlab("% epithelial cells")

![](HLCA_core_atlas_files/figure-markdown_strict/unnamed-chunk-20-1.png)

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
In the second [example](./HLCA_core_atlas_supervised.Rmd) we propose a
supervised workflow using this annotation to guide both the metacell
identification and the batch correction.

We can save the results for comparison with the second example.

    saveRDS(combined.mc,"./HLCA_data/combined.mc.unsup.rds")
