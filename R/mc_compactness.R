#' Compute compactness of metacells.
#'
#' \code{mc_compactness} 
#' This function computes the compactness of metacells (relies on python functions adapted from https://github.com/dpeerlab/SEACells/blob/main/SEACells/evaluate.py). 
#' @param cell.membership A data frame containing at least one column assigning each single-cell to a metacell (named "membership" if no other label is provided in group.label)
#' and single-cell IDs as rownames. 
#' @param sc.obj A Seurat object containing the single-cell data from which the metacells were built.
#' @param sc.reduction (optional, default is "pca") A string indicating which low embedding from sc.obj should be used to compute compactness or
#' a data frame containing a pre-computed embedding of the single-cell data. 
#' @param group.label (optional, default is "membership") A string indicating the column name from *cell.membership* that should be used to compute 
#' the compactness metric. 
#' @param dims (optional, default is NULL) Vector indicating the embedding components that should be used to compute the compactness. 
#' If NULL, components 1 to max number of components in the embedded space will be used.
#' @return  A vector containing the compactness of each metacell.
#' @examples
#' mc_compactness(cell.membership = CD34_mc@misc$cell_membership$membership, sc.obj = CD34_sc, sc.reduction = "pca")
#' @export
#' 
# mc_compactness(sc.obj = MetacellAnalysisToolkit::CD34_sc, sc.reduction = "pca", cell.membership = MetacellAnalysisToolkit::CD34_mc@misc$cell_membership)
# CD34_mc@meta.data["compactness"] <- mc_compactness(sc.obj = MetacellAnalysisToolkit::CD34_sc,
#                                                   sc.reduction = "pca",
#                                                   cell.membership = MetacellAnalysisToolkit::CD34_mc@misc$cell_membership)
# head(CD34_mc@meta.data)
# CD34_mc@meta.data["compactness"] <- mc_compactness(sc.obj = MetacellAnalysisToolkit::CD34_sc,
#                                                   sc.reduction = "pca",
#                                                   cell.membership = data.frame(
#                                                     row.names = rownames(MetacellAnalysisToolkit::CD34_mc@misc$cell_membership),
#                                                     grouping = MetacellAnalysisToolkit::CD34_mc@misc$cell_membership$membership
#                                                   ),
#                                                   group.label = "grouping")
# 
# head(CD34_mc@meta.data)

# diffusion_comp <- get_dim_reduc(sc.obj = CD34_sc, n.components = 30)
# CD34_mc@meta.data["compactaness"] <- mc_compactness(cell.membership = CD34_mc@misc$cell_membership, 
#                                                     sc.obj = CD34_sc, 
#                                                     sc.reduction = diffusion_comp)

# add option to provide a MC seurat and if provided add the compactness to the MC_seurat
mc_compactness <- function(cell.membership, sc.obj, sc.reduction, group.label = "membership", dims=NULL) {
  if(is.null(sc.obj)){
    if(!(is.data.frame(sc.reduction) | is.matrix(sc.reduction)) ){
      stop("sc.obj is NULL, sc.reduction should be a dataframe (or matrix) containing the components (columns) of single-cell embedding")
    }  
  }
  
  if(assertthat::is.string(sc.reduction)){
    # if sc_reduction does not exist compute pca:
    if(is.null(sc.obj@reductions[[sc.reduction]])){
      message("Low dimensionnal embedding not found in sc.obj")
      message("Computing PCA ...")
      sc.obj <- Seurat::NormalizeData(sc.obj, normalization.method = "LogNormalize")
      sc.obj <- Seurat::FindVariableFeatures(sc.obj, nfeatures = 2000)
      sc.obj <- Seurat::ScaleData(sc.obj)
      sc.obj <- Seurat::RunPCA(sc.obj, verbose = F)
      sc.reduction <- Seurat::Embeddings(sc.obj,reduction = sc.reduction)
    }else{
      sc.reduction <- Seurat::Embeddings(sc.obj,reduction = sc.reduction)
    }  
  } else if(!(is.data.frame(sc.reduction) | is.matrix(sc.reduction)) ){
    stop("sc.reduction should be a string indicating the name of the embedding to use in the reduction slot of sc.obj or a dataframe (or matrix) containing the components (columns) of single-cell embedding")
  }  
  
  if (is.null(dims)) {
    dims <- c(1:dim(sc.reduction)[2])
  }
  
  # remove MetaCell2 outliers 
  memberships_without_outliers <- na.exclude(cell.membership) 
  membership_vector <- memberships_without_outliers[, group.label]
  names(membership_vector) <- rownames(memberships_without_outliers)
  
  sc.reduction = sc.reduction[names(membership_vector), dims]
  
  centroids <- stats::aggregate(sc.reduction,
                                by = list(metacell = membership_vector),
                                FUN = var)
  compactness <- apply(centroids[, -1], 1, mean)
  names(compactness) <- centroids[, 1]
  
  return(compactness)
}
