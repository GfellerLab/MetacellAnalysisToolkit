#' Compute separation of metacells.
#'
#' \code{mc_separation} 
#' This function computes the separation of metacells (relies on python functions adapted from https://github.com/dpeerlab/SEACells/blob/main/SEACells/evaluate.py). 
#' @param cell.membership A data frame containing at least one column assigning each single-cell to a metacell (named "membership" if no other label is provided in group.label)
#' and single-cell IDs as rownames. 
#' @param sc.obj A Seurat object containing the single-cell data from which the metacells were built.
#' @param sc.reduction (optional, default is "pca") A string indicating which low embedding from sc.obj should be used to compute compactness or
#' a data frame containing a pre-computed embedding of the single-cell data. 
#' @param group.label (optional, default is "membership") A string indicating the column name from *cell.membership* that should be used to compute 
#' the compactness metric. 
#' @param dims (optional, default is NULL) Vector indicating the embedding components that should be used to compute the compactness. 
#' If NULL, components 1 to max number of components in the embedded space will be used.
#' @param nth.nbr (optional, default is 1) Number of nearest neighbors used to compute the separation metric. 
#' @return  A vector containing the separation metric for each metacell.
#' @examples
#' mc_separation(cell.membership = CD34_mc@misc$cell_membership, sc.obj = CD34_sc)
#' @export
#' 
# mc_separation(sc.obj = CD34_sc, sc.reduction = "pca", cell.membership = CD34_mc@misc$cell_membership)
# CD34_mc@meta.data["separation"] <- mc_separation(sc.obj = CD34_sc,
#                                                  sc.reduction = CD34_sc@reductions$pca@cell.embeddings, cell.membership = MetacellAnalysisToolkit::CD34_mc@misc$cell_membership)
# head(CD34_mc@meta.data)


# diffusion_comp <- get_dim_reduc(sc.obj = CD34_sc, n.components = 30)
# CD34_mc@meta.data["separation"] <- mc_separation(cell.membership = CD34_mc@misc$cell_membership,
#                                                    sc.obj = CD34_sc,
#                                                    sc.reduction = diffusion_comp, nth.nbr = 1)

# CD34_mc@meta.data["separation"] <- mc_separation(
#   sc.obj = CD34_sc,
#   sc.reduction = diffusion_comp,
#   cell.membership = data.frame(
#     row.names = rownames(CD34_mc@misc$cell_membership),
#     grouping = CD34_mc@misc$cell_membership$membership
#   ),
#   group.label = "grouping")


mc_separation <- function(cell.membership, sc.obj = NULL, sc.reduction = "pca", group.label = "membership", dims=NULL, nth.nbr = 1) {
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
                                FUN = mean)
  dist_matrix <- as.matrix(dist(centroids[-1]))
  separation_distances <- apply(dist_matrix, 1, function(x) {
    x[order(x)[nth.nbr + 1]]  
  })
  
  return(separation_distances)
}
