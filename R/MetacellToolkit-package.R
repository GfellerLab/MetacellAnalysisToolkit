#' CD34 dataset at the metacell level
#'
#' Metacells derived from a single-cell dataset of CD34 cells obtained from \href{https://www.nature.com/articles/s41587-023-01716-9}{Persad et al. (2023)} and available 
#' at the single-cell level in MetacellAnalysisToolkit::CD34_sc.
#' Metacells data were generated using the command lines based on SuperCell and described in \href{https://github.com/GfellerLab/MetacellAnalysisToolkit}{the MetacellAnalysisToolkit github repository}.
#'
#' @format A Seurat object containing the metacells data, *i.e.* gene expression matrix as well as the membership of each single-cell to a metacell in the *misc* slot of the Seurat object.
#' @source \url{https://www.nature.com/articles/s41587-023-01716-9}

"CD34_mc"


#' CD34 dataset at the single-cell level
#'
#' Single-cell dataset composed of 6881 CD34 cells obtained from \href{https://www.nature.com/articles/s41587-023-01716-9}{Persad et al. (2023)}.
#'
#' @format An Seurat object containing the single-cell data (raw counts and metadata).
#' @source \url{https://www.nature.com/articles/s41587-023-01716-9}

"CD34_sc"
