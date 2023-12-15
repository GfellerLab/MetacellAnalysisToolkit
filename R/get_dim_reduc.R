#' Compute separation of metacells.
#'
#' \code{get_dim_reduc} 
#' This function computes diffusion components from palantir starting from a Seurat object (relies on python functions adapted from https://github.com/dpeerlab/SEACells/blob/main/SEACells/evaluate.py).
#' @param sc.obj A Seurat object containing the single-cell data from which the metacells were built.
#' @param sc.reduction (optional, default is "pca") A string indicating which low embedding from sc.obj should be used to compute the diffusion 
#' components or a data frame containing a pre-computed embedding of the single-cell data.
#' @param n.components (optional, default is 10) Number of embedding components that should be used to compute diffusion maps. 
#' @return  A data.frame containing the diffusion components.
#' @examples
#' get_dim_reduc(sc.obj = CD34_sc, n.components = 30)
#' @export
#' 

# add option to provide a MC seurat and if provided add the diff components to the MC_seurat
get_dim_reduc <- function(sc.obj, sc.reduction = "pca", n.components = 10){
  
  reticulate::source_python(system.file("python/QC_functions.py", package = "MetacellAnalysisToolkit"))
  # reticulate::source_python("inst/python/QC_functions.py")
  
  if(assertthat::is.string(sc.reduction)){
    # if sc_reduction does not exist compute pca:
    if(is.null(sc.obj@reductions[[sc.reduction]])){
      message("Low dimensionnal embessing not found in sc.obj")
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
  
  
  # create single-cell anndata to compute separation using the python function
  sc_ad <- anndata::AnnData(
    X = Matrix::t(Seurat::GetAssayData(sc.obj, layer = "counts")),
    obsm = list(sc_reduction =  sc.reduction)
  )
  
  # compute separation 
  message("Computing diffusion maps ...")
  
  emb <- get_diffusion_map(
    sc_ad,
    low_dim_embedding = "sc_reduction",
    n_comp = as.integer(n.components)
  )
  
  return(emb)
} 

