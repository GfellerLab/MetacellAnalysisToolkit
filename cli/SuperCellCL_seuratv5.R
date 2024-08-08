library(Seurat)
library(getopt)
#library(doParallel)
# library(SuperCell)
if(packageVersion("Seurat") >= 5) {options(Seurat.object.assay.version = "v4"); print("you are using seurat v5 with assay option v4"
										
MetacellRawExpression_Seurat_v5 <- function (object, pb.method = "aggregate", assays = NULL, features = NULL, 
                                             return.seurat = TRUE, group.by = "ident", add.ident = NULL, 
                                             layer = "counts", verbose = TRUE, ...) 
{
  CheckDots(..., fxns = "CreateSeuratObject")
  if (!is.null(x = add.ident)) {
    .Deprecated(msg = "'add.ident' is a deprecated argument, please use the 'group.by' argument instead")
    group.by <- c("ident", add.ident)
  }
  if (!(pb.method %in% c("average", "aggregate"))) {
    stop("'pb.method' must be either 'average' or 'aggregate'")
  }
  object.assays <- FilterObjects(object = object, classes.keep = "Assay")
  assays <- assays %||% object.assays
  if (!all(assays %in% object.assays)) {
    assays <- assays[assays %in% object.assays]
    if (length(x = assays) == 0) {
      stop("None of the requested assays are present in the object")
    }
    else {
      warning("Requested assays that do not exist in object. Proceeding with existing assays only.")
    }
  }
  if (length(x = layer) == 1) {
    layer <- rep_len(x = layer, length.out = length(x = assays))
  }
  else if (length(x = layer) != length(x = assays)) {
    stop("Number of layers provided does not match number of assays")
  }
  data <- FetchData(object = object, vars = rev(x = group.by))
  data <- data[which(rowSums(x = is.na(x = data)) == 0), , 
               drop = F]
  if (nrow(x = data) < ncol(x = object)) {
    message("Removing cells with NA for 1 or more grouping variables")
    object <- subset(x = object, cells = rownames(x = data))
  }
  for (i in 1:ncol(x = data)) {
    data[, i] <- as.factor(x = data[, i])
  }
  num.levels <- sapply(X = 1:ncol(x = data), FUN = function(i) {
    length(x = levels(x = data[, i]))
  })
  if (any(num.levels == 1)) {
    message(paste0("The following grouping variables have 1 value and will be ignored: ", 
                   paste0(colnames(x = data)[which(num.levels <= 1)], 
                          collapse = ", ")))
    group.by <- colnames(x = data)[which(num.levels > 1)]
    data <- data[, which(num.levels > 1), drop = F]
  }
  if (ncol(x = data) == 0) {
    message("All grouping variables have 1 value only. Computing across all cells.")
    category.matrix <- matrix(data = 1, nrow = ncol(x = object), 
                              dimnames = list(Cells(x = object), "all"))
    if (pb.method == "average") {
      category.matrix <- category.matrix/sum(category.matrix)
    }
  }
  else {
    category.matrix <- Matrix::sparse.model.matrix(object = as.formula(object = paste0("~0+", 
                                                                                       paste0("data[,", 1:length(x = group.by), "]", collapse = ":"))))
    colsums <- colSums(x = category.matrix)
    category.matrix <- category.matrix[, colsums > 0]
    colsums <- colsums[colsums > 0]
    if (pb.method == "average") {
      category.matrix <- Sweep(x = category.matrix, MARGIN = 2, 
                               STATS = colsums, FUN = "/")
    }
    colnames(x = category.matrix) <- sapply(X = colnames(x = category.matrix), 
                                            FUN = function(name) {
                                              name <- gsub(pattern = "data\\[, [1-9]*\\]", 
                                                           replacement = "", x = name)
                                              return(paste0(rev(x = unlist(x = strsplit(x = name, 
                                                                                        split = ":"))), collapse = "_"))
                                            })
  }
  data.return <- list()
  for (i in 1:length(x = assays)) {
    data.use <- GetAssayData(object = object, assay = assays[i], 
                             layer = layer[i])
    features.to.avg <- features %||% rownames(x = data.use)
    if (inherits(x = features, what = "list")) {
      features.to.avg <- features[i]
    }
    if (IsMatrixEmpty(x = data.use)) {
      warning("The ", layer[i], " layer for the ", assays[i], 
              " assay is empty. Skipping assay.", immediate. = TRUE, 
              call. = FALSE)
      next
    }
    bad.features <- setdiff(x = features.to.avg, y = rownames(x = data.use))
    if (length(x = bad.features) > 0) {
      warning("The following ", length(x = bad.features), 
              " features were not found in the ", assays[i], 
              " assay: ", paste(bad.features, collapse = ", "), 
              call. = FALSE, immediate. = TRUE)
    }
    features.assay <- intersect(x = features.to.avg, y = rownames(x = data.use))
    if (length(x = features.assay) > 0) {
      data.use <- data.use[features.assay, ]
    }
    else {
      warning("None of the features specified were found in the ", 
              assays[i], " assay.", call. = FALSE, immediate. = TRUE)
      next
    }
    
    data.return[[i]] <- as.sparse(x = (data.use %*% category.matrix))
    colnames(data.return[[i]]) <- paste0("Metacell_",c(1:ncol(data.return[[i]])))
    
    names(x = data.return)[i] <- assays[[i]]
  }
  if (return.seurat) {
    toRet <- CreateSeuratObject(counts = data.return[[1]], 
                                project = if (pb.method == "average") 
                                  "Average"
                                else "Aggregate", assay = names(x = data.return)[1], 
                                ...)
    # toRet <- SetAssayData(object = toRet, assay = names(x = data.return)[1], 
    #                       layer = "data", new.data = log1p(x = as.matrix(x = data.return[[1]])))
    
    if (length(x = data.return) > 1) {
      for (i in 2:length(x = data.return)) {
        
        toRet[[names(x = data.return)[i]]] <- CreateAssay5Object(counts = data.return[[i]])
        # toRet <- SetAssayData(object = toRet, assay = names(x = data.return)[i], 
        #                       layer = "data", new.data = log1p(x = as.matrix(x = data.return[[i]])))
        
      }
    }
    if (DefaultAssay(object = object) %in% names(x = data.return)) {
      DefaultAssay(object = toRet) <- DefaultAssay(object = object)
    }
    if ("ident" %in% group.by) {
      first.cells <- c()
      for (i in 1:ncol(x = category.matrix)) {
        first.cells <- c(first.cells, Position(x = category.matrix[, 
                                                                   i], f = function(x) {
                                                                     x > 0
                                                                   }))
      }
      Idents(object = toRet) <- Idents(object = object)[first.cells]
    }
    return(toRet)
  }
  else {
    return(data.return)
  }
}      
)}

spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'input',  'i', 1, "character", "REQUIRED : Path to seurat object (.rds) or anndata (.h5ad) object",
  'outdir',     'o', 1, "character", 'Outdir path (default ./)',
  'isNorm', 'd', 0, "logical", "Whether the data are already normalised (default FALSE and log normalize data using Seurat)",
  "nPCs", "n", 1, "numeric", "Number of principal components to use for construction of single-cell kNN network (default 50)",
  "nFeatures", "f", 1, "numeric", "number of highly variable features to use (default 2000), selected with Seurat",
  "k.knn", "k", 1, "numeric", "k for the knn used in the wnn used for metacell identification (default 30)",
  "gamma", "g", 1, "numeric", "gamma used for metacell identification (default 10)",
  "annotations", "a", 1, "character", "Identify metacell per label regarding a gigven annotation (column name in the metadata of the object)",
  "cores", "l", 1, "numeric", "number of cores to use for parallel identifications of metacells in each label",
  "minMetacells", "m", 1, "numeric", "minimum number of metacells to identify (usefull if an annotation is given)",
  "output", "s", 1, "character", "how to save the SuperCell results ('seurat' : Seurat object (.rds), 'adata' : anndata object (.h5ad), 'SC' : SuperCell object (SC.rds))"
), byrow=TRUE, ncol=5)

opt = getopt(spec)

if ( !is.null(opt$help) | (is.null(opt$input) & is.null(opt$output))) {
  cat("Command line implemention of SuperCell workflow:")
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}



if (is.null(opt$outdir)) {
  opt$outdir <- "./"
}

if (is.null(opt$nPCs)) {
  opt$nPCs <- 50
}

if (is.null(opt$k.knn)) {
  opt$k.knn <- 30
}

if (is.null(opt$nFeatures)) {
  opt$nFeatures <- 2000
}

if (is.null(opt$gamma)) {
  opt$gamma <- 10
}

if (is.null(opt$isNorm)) {
  opt$isNorm <- FALSE
}


print(opt)

if(!is.null(opt$cores) & !is.null(opt$annotations)){
  
  library(doParallel)
  library(parallel)
  
  cluster <- parallel::makeCluster(opt$cores)
  
  doParallel::registerDoParallel(cluster)
}


dir.create(opt$outdir,recursive = T,showWarnings = F)

if (endsWith(x = opt$input,suffix = ".h5ad")) {
  
  adata <- anndata::read_h5ad(opt$input)
  if (!is.null(adata$raw$X)) {
    countMatrix <- Matrix::t(adata$raw$X)
  } else {
    countMatrix <-  Matrix::t(adata$X)
  }
  colnames(countMatrix) <- adata$obs_names
  rownames(countMatrix) <- adata$var_names
  if (!opt$isNorm) {
    sobj <- CreateSeuratObject(counts = countMatrix,meta.data = adata$obs)
    remove(countMatrix)
  } else { # we assume norm counts are in .X and raw counts in .raw.X
    normMatrix <-  Matrix::t(adata$X)
    colnames(normMatrix) <- adata$obs_names
    rownames(normMatrix) <- adata$var_names
    sobj <- CreateSeuratObject(counts = normMatrix,meta.data = adata$obs)
    sobj@assays$RNA@data <- sobj@assays$RNA@counts # countMatrix will be use for aggregation after metacell identification
    remove(normMatrix)
  }
  remove(adata)
  gc(verbose = F)
  
} else {
  sobj <- readRDS(opt$input)
}




if(is.null(opt$cores)){
  opt$cores = 1
} 

if (!opt$isNorm) {
  cat("Normalize data...")
  sobj <- NormalizeData(sobj,verbose = F)
}

fields <- sapply(X = colnames(sobj@meta.data) , 
                 FUN = function(X) {is.character(sobj[[X]][,1]) | is.factor(sobj[[X]][,1])})
fields <- names(fields[which(fields)])

cat("Identify Metacells...\n")


if (opt$cores > 1 & !is.null(opt$annotations)) {
  
  SCs <- foreach::foreach(sobj.label = SplitObject(sobj,split.by = opt$annotations)) %dopar% {
    
    #print(paste0("Treat ",label, " cells"))
    #sobj.label <- sobj[,sobj[[opt$annotations]][,1] == label]
    
    # adapt parameters regarding number of single cells (occurs mainly when an annotation is given)
    if (is.null(opt$minMetacells)) {
      opt$minMetacells <- 1
    }
    minMetacells <- min(ncol(sobj.label),opt$minMetacells)
    k.knn <- min(opt$k.knn, ncol(sobj.label)-1)
    n.pc <- min(opt$nPCs, ncol(sobj.label)-1)
    
    
    
    targetGamma <- min(ncol(sobj.label)/minMetacells,opt$gamma)
    
    
    if (ncol(sobj.label)>4) { 
    
    sobj.label <- FindVariableFeatures(sobj.label,nfeatures = opt$nFeatures,verbose = F) #is performed on raw counts (as in Seurat workflow) only if is.norm = F 
    
    
    
    cat(paste0("Identify ",round(ncol(sobj.label)/targetGamma)," metacells using SuperCell...\n"))
    
    SC.label <- SuperCell::SCimplify(GetAssayData(sobj.label,slot = "data"),  # normalized gene expression matrix 
                                     n.pc = n.pc,
                                     k.knn = k.knn, # number of nearest neighbors to build kNN network
                                     gamma = targetGamma, # graining level
                                     genes.use = VariableFeatures(sobj.label))
  } else {
    cat("object contain less than 5 single cells, simplification is not possible\naggregating all single-cells in one metacell.")
    membership <- c(rep(1,ncol(sobj.label)))
    names(membership) <- colnames(sobj.label)
    SC.label <- list("N.SC" = 1,
                     "membership" = membership,
                     "supercell_size" = c(ncol(sobj.label)))
  }
    
    #SCs[[label]] <- SC.label
  }
  
  parallel::stopCluster(cluster)
  gc(verbose = F)
  
} else {
  

  
  SCs <- list()
  if (is.null(opt$annotations)) {
    opt$annotations <- "orig.ident"
  }
  
  if(length(unique(sobj[[opt$annotations]][,1]))> 1) {
    cat("Identify Metacells sequentially...\n")
  }
  
  for (label in unique(sobj[[opt$annotations]][,1])) {
    
    if(length(unique(sobj[[opt$annotations]][,1]))> 1) {
      cat(paste0("Treat ",label, " cells\n"))
    }
    
    sobj.label <- sobj[,sobj[[opt$annotations]][,1] == label]
    
    # adapt parameters regarding number of single cells (occurs mainly when an annotation is given)
    if (is.null(opt$minMetacells)) {
      opt$minMetacells <- 1
    }
    minMetacells <- min(ncol(sobj.label),opt$minMetacells)
    k.knn <- min(opt$k.knn, ncol(sobj.label)-1)
    n.pc <- min(opt$nPCs, ncol(sobj.label)-1)
    
    
    
    targetGamma <- min(ncol(sobj.label)/minMetacells,opt$gamma)
    
    if (ncol(sobj.label)>4) { 
    
    sobj.label <- FindVariableFeatures(sobj.label,nfeatures = opt$nFeatures,verbose = F) #is performed on raw counts in Seurat
    
    
    fields <- sapply(X = colnames(sobj.label@meta.data) , 
                     FUN = function(X) {is.character(sobj.label[[X]][,1]) | is.factor(sobj.label[[X]][,1])})
    fields <- names(fields[which(fields)])
    
    cat(paste0("Identify ",round(ncol(sobj.label)/targetGamma)," metacells using SuperCell...\n"))
    if (!exists("normMatrix")) {
      SC.label <- SuperCell::SCimplify(GetAssayData(sobj.label,slot = "data"),  # normalized gene expression matrix 
                                       n.pc = n.pc,
                                       k.knn = k.knn, # number of nearest neighbors to build kNN network
                                       gamma = targetGamma, # graining level
                                       genes.use = VariableFeatures(sobj.label))
    } else {
      SC.label <- SuperCell::SCimplify(normMatrix,  # normalized gene expression matrix in case of normalized anndata object as input
                                       n.pc = n.pc,
                                       k.knn = k.knn, # number of nearest neighbors to build kNN network
                                       gamma = targetGamma, # graining level
                                       genes.use = VariableFeatures(sobj.label))
      remove(normMatrix)
      gc(verbose = F)
    }
    } else {
      cat("object contain less than 5 single cells, simplification is not possible\naggregating all single-cells in one metacell.")
      membership <- c(rep(1,ncol(sobj.label)))
      names(membership) <- colnames(sobj.label)
      SC.label <- list("N.SC" = 1,
                       "membership" = membership,
                       "supercell_size" = c(ncol(sobj.label)))
    }
    
    SCs[[label]] <- SC.label
  }
}

SC <- SuperCell::supercell_merge(SCs)

if (opt$output != "SC") {
  
  if (opt$isNorm) {
    # in case normalized data was used and assigned to counts and data slot prior metacell identification,
    # we need to return to original raw counts for data aggregation
    sobj <- CreateSeuratObject(counts = countMatrix,meta.data = sobj@meta.data) 
    gc(verbose = F)
  }
  
  print(options("Seurat.object.assay.version"))
  
  sobj$Metacell <- SC$membership
  if (packageVersion("Seurat") < 5) {
  sobjMC <- AggregateExpression(sobj,
                                        return.seurat = T,
                                        group.by = "Metacell",
                                        slot = "counts", 
                                        verbose  = F)
  sobjMC@assays$RNA@data <- sobjMC@assays$RNA@counts # because AggregateExpression with slot = counts set data to log1p(aggregated counts)
  } else {
      sobjMC <- MetacellRawExpression_Seurat_v5(seurat,
                                                group.by = "Metacell",
                                                return.seurat = T)
      
      sobjMC[["RNA"]] <- as(object = sobjMC[["RNA"]], Class = "Assay")
  }
  gc(verbose = F)
  
  
  cat("Assign metadata to metacells and compute purities...\n")
  for (f in fields) {
    
    if (!all(is.na(sobj[[f]][,1]))) {
      sobjMC[[f]] <- SuperCell::supercell_assign(clusters = sobj[[f]][SC$cell.ids,1],
                                                 supercell_membership =  SC$membership,
                                                 method = "absolute")
      if(length(unique(sobj[[f]][,1]))> 1) {
        sobjMC[[paste0(f,"_purity")]] <- SuperCell::supercell_purity(clusters = sobj[[f]][SC$cell.ids,1],
                                                                     supercell_membership =  SC$membership)
      }
    }
  }
  
  #Store size
  sobjMC$size <- SC$supercell_size
  
  #Store membership in sobjMC
  membershipDF <- data.frame(SC$membership)
  colnames(membershipDF) <- "membership"
  sobjMC@misc$cell_membership <- membershipDF
  
  
  if (opt$output == "adata") {
    adataMC <- anndata::AnnData(X = Matrix::t(GetAssayData(object = sobjMC,slot = "counts")),
                                obs = sobjMC@meta.data,
                                uns =  list(cell_membership = membershipDF))
    anndata::write_h5ad(anndata = adataMC,filename = paste0(opt$outdir,"/mc_adata.h5ad"))
  }
  
  if (opt$output == "seurat") {
    saveRDS(object = sobjMC, file = paste0(opt$outdir,"/mc_Seurat.rds"))
  }
  
} else {
  saveRDS(object = SC, file = paste0(opt$outdir,"/mc_SuperCell.rds"))
}
































