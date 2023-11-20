#' Compute compactness of metacells.
#'
#' \code{mc_compactness} 
#' This function computes the compactness of metacells (relies on python functions adapted from https://github.com/dpeerlab/SEACells/blob/main/SEACells/evaluate.py). 
#' @param cell.membership A data frame containing at least one column assigning each single-cell to a metacell (named "membership" if no other label is provided in group.label)
#' and single-cell IDs as rownames. 
#' @param sc.obj A Seurat object containing the single-cell data from which the metacells were built.
#' @param sc.reduction (optional, default is "pca") A string indicating which low embedding from sc.obj should be used to compute compactness or
#' a data frame containing a pre-computed embedding of the single-cell data. If sc.reduction is a string and sc.reduction is not found in the sc.obj embeddings
#' PCA will be performed on the single-cell data.
#' @param group.label (optional, default is "membership") A string indicating the column name from *cell.membership* that should be used to compute the compactness metric. 
#' @param diffusion.components (optional, default is TRUE) A boolean indicating if compactness should be computed on diffusion components computed based on the low 
#' dimensionnal embedding provided in *sc.reduction*.
#' @param n.components (optional, default is 30) Number of embedding components that should be used to compute the compactness. 
#' @param n.features (optional, default is 2000) Number of features that should be used to run PCA if sc.reduction is not found. 
#' @return  A vector containing the compactness of each metacell.
#' @examples
#' mc_compactness(cell.membership = CD34_mc@misc$cell_membership, sc.obj = CD34_sc)
#' @export
#' 
# mc_compactness(sc.obj = MetacellAnalysisToolkit::CD34_sc, sc.reduction = "pca", cell.membership = MetacellAnalysisToolkit::CD34_mc@misc$cell_membership)
# CD34_mc@meta.data["compactness"] <- mc_compactness(sc.obj = MetacellAnalysisToolkit::CD34_sc,
#                                                   sc.reduction = CD34_sc@reductions$pca@cell.embeddings,
#                                                   cell.membership = MetacellAnalysisToolkit::CD34_mc@misc$cell_membership)
# head(CD34_mc@meta.data)
# CD34_mc@meta.data["compactness"] <- mc_compactness(sc.obj = MetacellAnalysisToolkit::CD34_sc,
#                                                   sc.reduction = CD34_sc@reductions$pca@cell.embeddings,
#                                                   cell.membership = data.frame(
#                                                     row.names = rownames(MetacellAnalysisToolkit::CD34_mc@misc$cell_membership),
#                                                     grouping = MetacellAnalysisToolkit::CD34_mc@misc$cell_membership$membership
#                                                   ),
#                                                   group.label = "grouping")
# 
# head(CD34_mc@meta.data)

mc_compactness <- function(cell.membership, sc.obj, sc.reduction = "pca", group.label = "membership", diffusion.components = TRUE, n.components = 30, n.features = 2000){
  reticulate::source_python(system.file("python/QC_functions.py", package = "MetacellAnalysisToolkit"))

  if(assertthat::is.string(sc.reduction)){
    # if sc_reduction does not exist compute pca:
    if(is.null(sc.obj@reductions[[sc.reduction]])){
      message("Low dimensionnal embessing not found in sc.obj")
      message("Computing PCA ...")
      sc.obj <- Seurat::NormalizeData(sc.obj, normalization.method = "LogNormalize")
      sc.obj <- Seurat::FindVariableFeatures(sc.obj, nfeatures = n.features)
      sc.obj <- Seurat::ScaleData(sc.obj)
      sc.obj <- Seurat::RunPCA(sc.obj, verbose = F)
      sc.reduction <- Seurat::Embeddings(sc.obj@reductions[[sc.reduction]])
    }else{
      sc.reduction <- Seurat::Embeddings(sc.obj@reductions[[sc.reduction]])
    }  
  } else if(!(is.data.frame(sc.reduction) | is.matrix(sc.reduction)) ){
    stop("sc.reduction should be a string indicating the name of the embedding to use in the reduction slot of sc.obj or a dataframe (or matrix) containing the components (columns) of single-cell embedding")
  }  
  
  # Complete membership if not all single-cell found in the membership data frame
  cell.membership[colnames(sc.obj)[!colnames(sc.obj) %in% rownames(cell.membership)],] <- NA
  
  # create single-cell anndata to compute compactness using the python function
  sc_ad <- anndata::AnnData(
    X = Matrix::t(Seurat::GetAssayData(sc.obj, slot = "counts")),
    obs = cell.membership,
    obsm = list(sc_reduction =  sc.reduction)
    )
  
  # compute compactness 
  message("Computing compactness ...")
  metric_res <- compactness(
    sc_ad,
    low_dim_embedding = 'sc_reduction',
    MC_label = group.label,
    DO_DC = diffusion.components,
    name = 'compactness',
    n_comp = as.integer(n.components)
    )
  
  return(metric_res$compactness)
} 


