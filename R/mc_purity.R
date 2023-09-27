#' Compute purity of metacells.
#'
#' \code{mc_purity} 
#' This function computes the purity of metacells. (code has been retrieved from the SuperCell::supercell_purity function).
#' @param membership: Membership vector associating each single-cell to a metacell. 
#' @param annotation: Vector of single-cell annotations of length equal to the number of single-cells used 
#' to build the metacell object. The annotations order should match the single-cell membership provided in the membership parameter.
#' @param method: Method to compute metacell purity. "max_proportion" if the purity is defined as a proportion of the most abundant annotation group (e.g. cell type) 
#' within the metacell or "entropy" if the purity is defined as the Shanon entropy of the annotation groups metacells consists of. 
#' @return  A vector containing the purity of each metacell.
#' @examples
#' mc_purity(membership = CD34_mc@misc$cell_membership$membership, annotation = CD34_sc$celltype)
#' @export


mc_purity <- function (membership, annotation, method = c("max_proportion", "entropy")[1]) {
  if (!(method %in% c("max_proportion", "entropy"))) {
    stop(paste("Method", method, "is not known. The available methods are:", 
               paste(method, collapse = ",")))
  }
  cl.gr <- table(annotation, membership)
  switch(method, entropy = {
    res <- apply(cl.gr, 2, entropy::entropy)
  }, max_proportion = {
    cluster.size <- as.numeric(table(annotation))
    group.size <- as.numeric(table(membership))
    Ng <- length(group.size)
    group.max.cl <- rep(0, Ng)
    cl.gr <- sweep(cl.gr, 2, group.size, "/")
    res <- apply(cl.gr, 2, max)
  }, {
    stop(paste("Method", method, "is not known. The available methods are:", 
               paste(method, collapse = ",")))
  })
  return(res)
}