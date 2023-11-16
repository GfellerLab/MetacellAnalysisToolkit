#' mc_projection of metacells in the single cell space  
#' 
#' \code{mc_projection} 
#' This function plots metacells in a single-cell space (pca, umap,..) taking the average coordinates of single cells in each metacell. 
#' @param sc.obj A Seurat object containing the single-cell data from which the metacells were built.
#' @param mc.obj A Seurat object containing the metacells data. If NULL, memberhsip should be provided as a data frame.
#' @param cell.membership A data frame containing at least one column assigning each single-cell to a metacell (named "membership" if no other label is provided in group.label)
#' and single-cell IDs as rownames. 
#' If *cell.membership* is NULL, the membership information will be retrieved from the *misc* slot in the *mc.obj* under the assumption that the metacell Seurat object
#' was generated using the MetacellToolkit command lines (see format of the *MetacellToolkit::CD34_mc* object).
#' @param metacell.label String corresponding to the name of the metadata column in mc.obj that should be used to color metacells. 
#' If mc.obj is NULL, this parameter will be ignored.
#' @param sc.label String corresponding to the name of the metadata column in sc.obj that should be used to color single cells.
#' @param dims Numerical vector (of length 2) indicating which components should be used in the 2D visualization of the data.
#' @param sc.reduction (optional, default is "umap") Either a string indicating which low embedding from sc.obj should be used for the 2D visualization of the data
#' or a data frame containing a pre-computed embedding of the single-cell data. If *sc.reduction* is a string which does not correspond to any embedding from sc.obj
#' PCA will be performed on the single-cell data and UMAP will be run to obtain the 2D representation. 
#' @param mc.color Colors for metacell idents
#' @param sc.color Colors for single-cell idents
#' @param alpha Transparency value for the single-cell points
#' @param pt_size Size the single-cell points
#' @param metric Column name of a continuous metric to use either to color the metacells (if continuous_metric is TRUE) or to define metacells sizes (if continuous_metric is FALSE).
#' @param continuous_metric Bolean indicating if the metric variable is continuous or not. If TRUE a continuous color scale will be used for the metacell colors.
#' @return A plot of metacells projected in the single cell space.
#' @examples
#' mc_projection(sc.obj = CD34_sc, mc.obj = CD34_mc)
#' mc_projection(sc.obj = CD34_sc, mc.obj = CD34_mc, metric = "celltype_purity", continuous_metric = TRUE)
#' mc_projection(sc.obj = CD34_sc, mc.obj = CD34_mc, metric = "celltype_purity", continuous_metric = TRUE, sc.label = "celltype")
#' mc_projection(sc.obj = CD34_sc, mc.obj = CD34_mc, sc.label = "celltype", metacell.label = "celltype")
#' 
#' @export

# tests:
# CD34_sc <- Seurat::NormalizeData(CD34_sc, normalization.method = "LogNormalize")
# CD34_sc <- Seurat::FindVariableFeatures(CD34_sc, nfeatures = 1000)
# CD34_sc <- Seurat::ScaleData(CD34_sc)
# CD34_sc <- Seurat::RunPCA(CD34_sc)
# CD34_sc <- Seurat::RunUMAP(CD34_sc, reduction = "pca", dims = c(1:10), n.neighbors = 10)
# 
# mc_projection(sc.obj = CD34_sc, mc.obj = CD34_mc, sc.reduction = "umap")
# mc_projection(sc.obj = CD34_sc, mc.obj = CD34_mc, metric = "celltype_purity", continuous_metric = T)
# mc_projection(sc.obj = CD34_sc, mc.obj = CD34_mc, metric = "celltype_purity", continuous_metric = T, sc.label = "celltype")
# mc_projection(sc.obj = CD34_sc, mc.obj = CD34_mc, metacell.label = "celltype")
# mc_projection(sc.obj = CD34_sc, mc.obj = CD34_mc, sc.label = "celltype", metacell.label = "celltype")
# 
# colors <- c("HMP" = "#BC80BD", "DCPre" = "#66A61E", "HSC" = "#A6761D",
#             "Ery" = "#E41A1C", "MEP" = "#B3B3B3", "cDC" = "#A6D854", "CLP" = "#1F78B4", "Mono" = "#E6AB02", "pDC" = "#B2DF8A" )
# mc_projection(sc.obj = CD34_sc, mc.obj = CD34_mc, sc.label = "celltype", metacell.label = "celltype", sc.color = colors, mc.color = colors)
# mc_projection(sc.obj = CD34_sc, mc.obj = CD34_mc, sc.label = "celltype", metacell.label = "celltype", sc.color = colors)
# mc_projection(sc.obj = CD34_sc, mc.obj = CD34_mc, sc.label = "celltype", metacell.label = "celltype", mc.color = colors)

mc_projection <- function(sc.obj,
                          mc.obj = NULL,
                          cell.membership = NULL,
                          metacell.label = NULL,
                          sc.label = NULL,
                          dims = c(1, 2),
                          sc.reduction = "umap",
                          mc.color = NULL,
                          sc.color = NULL,
                          alpha = 1,
                          pt_size = 0,
                          metric = "size",
                          continuous_metric = F) {
  
  if(is.null(mc.obj) & is.null(cell.membership)){
    stop("A membership vector should be provided either through the membership parameter or should be available in mc.obj as described in the documentation.")
  } 
  if(is.null(cell.membership)){
    membership <- mc.obj@misc$cell_membership$membership
    names(membership) <- rownames(mc.obj@misc$cell_membership)
  } else {
    membership <- cell.membership$membership
    names(membership) <- rownames(cell.membership)
  }  
  
  sc.obj$Metacell <- membership

  if(assertthat::is.string(sc.reduction)){
    # if sc_reduction does not exist compute pca and run UMAP:
    if(is.null(sc.obj@reductions[[sc.reduction]])){
      message("Low dimensionnal embessing not found in sc.obj")
      message("Computing PCA ...")
      sc.obj <- Seurat::NormalizeData(sc.obj, normalization.method = "LogNormalize")
      sc.obj <- Seurat::FindVariableFeatures(sc.obj, nfeatures = 2000)
      sc.obj <- Seurat::ScaleData(sc.obj)
      sc.obj <- Seurat::RunPCA(sc.obj, verbose = F)
      message("Running UMAP ...")
      sc.obj <- Seurat::RunUMAP(sc.obj, reduction = "pca", dims = c(1:30), n.neighbors = 15, verbose = F)
      scCoord <- Seurat::Embeddings(sc.obj@reductions[["umap"]])
    } else{
      scCoord <- Seurat::Embeddings(sc.obj@reductions[[sc.reduction]])
    } 
  } else if(!(is.data.frame(sc.reduction) | is.matrix(sc.reduction)) ){
    stop("sc.reduction should be a string indicating the name of the embedding to use in the reduction slot of sc.obj or a dataframe (or matrix) containing the components (columns) of single-cell embedding")
  } else{
    scCoord <- sc.reduction
  }  
  
  scCoordMetacell <-  cbind(scCoord, membership)
  
  centroids <- stats::aggregate(scCoord~membership, scCoord, mean) #should be taken from object slot
  
  # if metacell.label and metric provided, mc.obj is mandatory if mc.obj not found membership mandatory and metric set to size which we compute from membership 
  # just add a message to warn the user 
  centroids[[metric]] <- mc.obj[[metric]][,1]
  
  if(is.null(metacell.label)) {
    metacell.label <- "MC"
    centroids[[metacell.label]] <- rep("red", length(mc.obj[[metric]]))
  } else {
    centroids[[metacell.label]] <- as.factor(mc.obj[[metacell.label]][,1])
    
    if(!is.null(sc.label)){
      if(metacell.label == sc.label){
        sc.obj@meta.data[, metacell.label] <- as.factor(sc.obj@meta.data[, metacell.label])
        centroids[[metacell.label]] <- factor(centroids[[metacell.label]], levels = levels(sc.obj@meta.data[, metacell.label]))
      } 
    } 
  }
  
  if(!is.null(sc.label)){
    scCoord <- data.frame(scCoord)
    scCoord[[sc.label]] <- sc.obj[[sc.label]][,1]
    
    p <- ggplot2::ggplot(scCoord,
                         ggplot2::aes_string(colnames(scCoord)[dims[1]],
                                             colnames(scCoord)[dims[2]],
                                             color = sc.label)) +
      ggplot2::geom_point(size=pt_size, alpha = alpha) +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 2)))
    
    
  }else{
    p <- ggplot2::ggplot(data.frame(scCoord),
                         ggplot2::aes_string(colnames(scCoord)[dims[1]],
                                             colnames(scCoord)[dims[2]])) +
      ggplot2::geom_point(size=pt_size, color = "grey", alpha = 1)
  } 
  
  if(!continuous_metric){
    p <- p + ggplot2::geom_point(data=centroids, ggplot2::aes_string(colnames(centroids)[1 + dims[1]],
                                                                     colnames(centroids)[1 + dims[2]],
                                                                     fill = metacell.label, size = metric), colour="black", pch=21) 
  }else{
    p <- p + ggplot2::geom_point(data=centroids, ggplot2::aes_string(colnames(centroids)[1 + dims[1]],
                                                                     colnames(centroids)[1 + dims[2]],
                                                                     fill = metric), colour="black", pch=21, size=2) 
  } 
  
  if(!is.null(metacell.label) & !is.null(sc.label) & !is.null(mc.color)){
    if(metacell.label == sc.label & !is.null(mc.color)){
      sc.color = mc.color
    } 
  } 
  
  if(!is.null(metacell.label) & !is.null(sc.label) & !is.null(sc.color)){
    if(metacell.label == sc.label & !is.null(sc.color)){
      mc.color = sc.color
    } 
  } 
  
  if (!is.null(mc.color) & !continuous_metric) {
    p <- p + ggplot2::scale_fill_manual(values = mc.color) +  ggplot2::theme_classic() 
  }
  if (!is.null(sc.color)) {
    p <- p + ggplot2::scale_color_manual(values = sc.color) +  ggplot2::theme_classic()
  }
  return(p)
}
