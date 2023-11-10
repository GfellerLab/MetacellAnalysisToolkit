#' Generate boxplot to visualize QC metrics distributions.
#'
#' \code{qc_boxplot} 
#' This function generate boxplot to visualize QC metrics distributions across metacells.
#' @param mc.obj A Seurat object containing the metacells data .
#' @param qc.metrics Vector of strings indicating which QC metric should be considered in the mc.obj metadata. 
#' @param split.by (optional): String indicating if the boxplot should be splitted based a metacell annotation available in the mc.obj metadata dataframe. 
#' @param y.lim (optional): Y axis limit. 
#' By default, the orig.ident variable is used.
#' @return Boxplots representing the distribution of the qc.metrics.
#' @examples
#' qc_boxplot(mc.obj = CD34_mc, qc.metrics = c("size", "celltype_purity"))
#' @export
#' @importFrom dplyr %>%
#' 

qc_boxplot <- function(mc.obj, qc.metrics = NULL, split.by = "orig.ident", y.lim = NULL){
  
  if(length(qc.metrics) == 1){
    if(split.by == "orig.ident"){
      bxp <- ggplot2::ggplot(mc.obj@meta.data, ggplot2::aes(x = get(split.by), y = get(qc.metrics))) +
        ggplot2::geom_boxplot() + ggplot2::ylab(qc.metrics) + ggplot2::xlab("Metacells") + ggplot2::theme(axis.text.x = ggplot2::element_blank())
    } else{
      bxp <- ggplot2::ggplot(mc.obj@meta.data, ggplot2::aes(x = get(split.by), y = get(qc.metrics))) +
        ggplot2::geom_boxplot() + ggplot2::ylab(qc.metrics) + ggplot2::xlab("Metacells") + 
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1))
    } 
    
    if (!is.null(y.lim)) {
      bxp <- bxp + ggplot2::ylim(y.lim)
    } 
  }else{
    library(dplyr)
    bxp <- mc.obj@meta.data %>%
      dplyr::select(qc.metrics[qc.metrics %in% colnames(mc.obj@meta.data)]) %>%
      tidyr::gather(na.rm = TRUE) %>%
      ggplot2::ggplot(ggplot2::aes(y = value, x = 0)) +
      ggplot2::geom_boxplot() +
      # geom_density(aes(x = value, y = stat(scaled)), inherit.aes = FALSE) +
      ggplot2::facet_wrap(~key, scales = 'free') + ggplot2::theme(axis.title.x = ggplot2::element_blank(), 
                                                                  axis.text.x = ggplot2::element_blank(), 
                                                                  axis.ticks = ggplot2::element_blank())
  } 
  
  return(bxp)
} 
