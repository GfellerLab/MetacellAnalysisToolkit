#' Compute inner-normalized variance (INV) of metacells.
#'
#' \code{mc_INV} 
#' This function computes the inner-normalized variance (INV) of metacells (relies on python functions inspired from https://github.com/tanaylab/metacells/blob/master/metacells/tools/quality.py).
#' @param sc.obj A Seurat object containing the single-cell data from which the metacells were built.
#' @param cell.membership A data frame containing at least one column assigning each single-cell to a metacell (or other cell grouping) and 
#' single-cell IDs as rownames.
#' @param group.label (optional, default is "membership") A string indicating the column name from *cell.membership* that should be used to compute the INV metric. 
#' @param assay (optional, default is "RNA") A string indicating the assay from *sc.obj* that should be used to compute the INV metric. 
#' Waning: if counts and data slots are identical, the counts data will be normalized using Seurat ("LogNormalize" method) and the normalized data will 
#' be used to compute INV values.
#' @return  A vector containing the INV value of each metacell.
#' @examples
#' mc_INV(sc.obj = MetacellToolkit::CD34_sc, cell.membership = CD34_mc@misc$cell_membership)
#' @export

mc_INV <- function(sc.obj, 
                   cell.membership, 
                   assay = "RNA",
                   group.label = "membership"){
  reticulate::source_python(system.file("python/QC_functions.py", package = "MetacellToolkit"))
  
  if(identical(sc.obj[[assay]]@counts@x, sc.obj[[assay]]@data@x)){
    message("Counts and data slots are identical.")
    message("Normalizing data ...")
    sc.obj <- Seurat::NormalizeData(sc.obj, normalization.method = "LogNormalize")
  } 
  
  # create single-cell anndata to compute INV using the python function
  sc_ad <- anndata::AnnData(
    X = Matrix::t(Seurat::GetAssayData(sc.obj, slot = "data")),
    obs = cell.membership
  )
  # sc_ad$obs[group.label] <- cell.membership[, group.label]

  # compute separation
  message("Computing INV ...")
  inv_val <- mc_inner_normalized_var(
    sc_ad,
    MC_label = group.label
  )

  # save INV vector in the metacell Seurat object
  # mc.obj@meta.data["INV"] <- metric_res$INV
  # library(dplyr)
  # mean_vals <- as.matrix(Seurat::GetAssayData(sc.obj, slot = "data")) %>%
  #   group_by(cell.membership$membership) %>%
  #   summarise_each(funs(mean(., na.rm = TRUE)))
  # Matrix::rowVars(sparse_mat)
  # gene_values <- 
  # inv_val <- stats::quantile(gene_values, 0.95)
  return(inv_val$`0.95`)
} 
