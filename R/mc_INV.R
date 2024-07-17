#' Compute inner-normalized variance (INV) of metacells.
#'
#' \code{mc_INV} 
#' This function computes the inner-normalized variance (INV) of metacells (relies on functions inspired from https://github.com/tanaylab/metacells/blob/master/metacells/tools/quality.py).
#' @param sc.obj A Seurat object containing the single-cell data from which the metacells were built.
#' @param cell.membership A data frame containing at least one column assigning each single-cell to a metacell (named "membership" if no other label is provided in group.label)
#' and single-cell IDs as rownames.
#' @param group.label (optional, default is "membership") A string indicating the column name from *cell.membership* that should be used to compute the INV metric. 
#' @param assay (optional, default is "RNA") A string indicating the assay from *sc.obj* that should be used to compute the INV metric. 
#' @param slot (optional, default is "counts") A string indicating the slot/layer from *sc.obj* that should be used to compute the INV metric. 
#' @param do.norm (optional, default is T) Logical value indicating whether the input slot data should be normalized (normalization performed within each metacell). 
#' @param NA_val (optional, default is 1) Numeric value indicating which value should be used to replace NA values. 
#' @param scaling_factor (optional, default is NULL) Used only if *do.norm* is TRUE. Numeric value indicating which scaling factor should be used in the counts normalization. 
#' If NULL, the median of the total number of UMIs in the cells belonging to the metacell is used. 
#' @return  A vector containing the INV value of each metacell.
#' @examples
#' mc_INV(sc.obj = MetacellAnalysisToolkit::CD34_sc, cell.membership = MetacellAnalysisToolkit::CD34_mc@misc$cell_membership)
#' @export

mc_INV <- function(sc.obj,
                   cell.membership,
                   group.label = "membership",
                   assay = "RNA",
                   slot = "counts",
                   do.norm = T,
                   NA_val = 1,
                   scaling_factor = NULL){
  
  # compute separation
  message("Computing INV ...")
  
  # remove MetaCell2 outliers 
  memberships_without_outliers <- na.exclude(cell.membership) 
  membership_vector <- memberships_without_outliers[, group.label]
  names(membership_vector) <- rownames(memberships_without_outliers)
  
  sc.obj <- sc.obj[,names(membership_vector)]
  
  get_INV_val <- function(x, norm) {
    if(norm){
      median_val <- median(Matrix::rowSums(x))
      x_normalized <- (x / Matrix::rowSums(x)) * ifelse(is.null(scaling_factor), median_val, scaling_factor)
      result <- (1 / Matrix::colMeans(x_normalized)) * sparseMatrixStats::colVars(x_normalized)
    }else{
      result <- (1 / Matrix::colMeans(x)) * sparseMatrixStats::colVars(x)
    }
    return(result)
  }
  
  data_matrix <- Matrix::t(Seurat::GetAssayData(sc.obj, assay = assay, slot = slot))

  result_list <- tapply(1:nrow(data_matrix), membership_vector, function(indices) {
    x <- data_matrix[indices, ]
    get_INV_val(x, norm = do.norm)
  })
  INV_val <- do.call(rbind, result_list)
  
  INV_val[is.na(INV_val)] <- NA_val
  
  # Compute the 95th percentile (quantile)
  INV_val_qt <- apply(INV_val, 1, function(x) quantile(x, 0.95, na.rm = TRUE))
  
  return(INV_val_qt)
}


