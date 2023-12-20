#' Compute inner-normalized variance (INV) of metacells.
#'
#' \code{mc_INV} 
#' This function computes the inner-normalized variance (INV) of metacells (relies on python functions inspired from https://github.com/tanaylab/metacells/blob/master/metacells/tools/quality.py).
#' @param sc.obj A Seurat object containing the single-cell data from which the metacells were built.
#' @param cell.membership A data frame containing at least one column assigning each single-cell to a metacell (named "membership" if no other label is provided in group.label)
#' and single-cell IDs as rownames.
#' @param group.label (optional, default is "membership") A string indicating the column name from *cell.membership* that should be used to compute the INV metric. 
#' @param assay (optional, default is "RNA") A string indicating the assay from *sc.obj* that should be used to compute the INV metric. 
#' Waning: if counts and data slots are identical, the counts data will be normalized using Seurat ("LogNormalize" method) and the normalized data will 
#' be used to compute INV values.
#' @return  A vector containing the INV value of each metacell.
#' @examples
#' mc_INV(sc.obj = MetacellAnalysisToolkit::CD34_sc, cell.membership = MetacellAnalysisToolkit::CD34_mc@misc$cell_membership)
#' @export

mc_INV <- function(sc.obj, 
                   cell.membership, 
                   assay = "RNA",
                   group.label = "membership"){
  if(identical(sc.obj[[assay]]@counts@x, sc.obj[[assay]]@data@x)){
    message("Counts and data slots are identical.")
    message("Normalizing data ...")
    sc.obj <- Seurat::NormalizeData(sc.obj, normalization.method = "LogNormalize")
  } 

  # compute separation
  message("Computing INV ...")
  # remove MetaCell2 outliers 
  memberships_without_outliers <- na.exclude(cell.membership) 
  membership_vector <- memberships_without_outliers[, group.label]
  names(membership_vector) <- rownames(memberships_without_outliers)
  
  sc.obj <- sc.obj[,names(membership_vector)]

  INV_val <- stats::aggregate(Matrix::t(Seurat::GetAssayData(sc.obj, slot = "data")),
                              by = list(metacell = membership_vector),
                              FUN = function(x) (1 / mean(x)) * var(x))
  rownames(INV_val) <- INV_val[, 1]
  INV_val <- INV_val[,-1]
  
  INV_val[is.na(INV_val)] <- 0
  
  # Compute the 95th percentile (quantile)
  INV_val_qt <- apply(INV_val, 1, function(x) quantile(x, 0.95, na.rm = TRUE))
  
  return(INV_val_qt)
} 

