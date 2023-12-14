install.packages("Seurat") # v5 not yet on conda
remotes::install_github("GfellerLab/SuperCell",upgrade = "never")
remotes::install_github("GfellerLab/MetacellAnalysisToolkit",upgrade = "never")
remotes::install_github("rstudio/reticulate",upgrade = "never")  #temporary fix for reading sparse matrix with R anndata https://github.com/rstudio/reticulate/issues/1417
