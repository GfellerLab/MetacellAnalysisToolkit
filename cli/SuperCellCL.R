library(Seurat)
library(getopt)
#library(doParallel)
# library(SuperCell)
if(packageVersion("Seurat") >= 5) {options(Seurat.object.assay.version = "v3"); message("you are using seurat v5 with assay option v3")}

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

if(is.null(opt$cores)){
  opt$cores = 1
} 


print(opt)




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
    if(packageVersion("Seurat") >= 5) { sobj[["RNA"]] <- as(object = sobj[["RNA"]], Class = "Assay") }
    remove(countMatrix)
  } else { # we assume norm counts are in .X and raw counts in .raw.X
    normMatrix <-  Matrix::t(adata$X)
    colnames(normMatrix) <- adata$obs_names
    rownames(normMatrix) <- adata$var_names
    sobj <- CreateSeuratObject(counts = normMatrix,meta.data = adata$obs)
    if(packageVersion("Seurat") >= 5) { sobj[["RNA"]] <- as(object = sobj[["RNA"]], Class = "Assay") }
    sobj@assays$RNA@data <- sobj@assays$RNA@counts # countMatrix will be use for aggregation after metacell identification
    remove(normMatrix)
  }
  remove(adata)
  gc(verbose = F)
  
} else {
  sobj <- readRDS(opt$input)
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
  library(doParallel)
  library(parallel)
  
  cluster <- parallel::makeCluster(opt$cores)
  
  doParallel::registerDoParallel(cluster)
  print('do for each')
  SCs <- foreach::foreach(sobj.label = SplitObject(sobj,split.by = opt$annotations)) %dopar% {
    
    
    # adapt parameters regarding number of single cells (occurs mainly when an annotation is given)
    if (is.null(opt$minMetacells)) {
      opt$minMetacells <- 1
    }
    minMetacells <- min(ncol(sobj.label),opt$minMetacells)
    k.knn <- min(opt$k.knn, ncol(sobj.label)-1)
    n.pc <- min(opt$nPCs, ncol(sobj.label)-1)
    
    
    
    targetGamma <- min(ncol(sobj.label)/minMetacells,opt$gamma)
    
    
    if (ncol(sobj.label)>4) { 
      
      sobj.label <- Seurat::FindVariableFeatures(sobj.label,nfeatures = opt$nFeatures,verbose = F) #is performed on raw counts (as in Seurat workflow) only if is.norm = F 
      if(packageVersion("Seurat") >= 5) {
        genes_exclude <- rownames(sobj.label)[which(Matrix::rowSums(Seurat::GetAssayData(sobj.label, layer = "data") != 0) < 2)] 
      } else {
        genes_exclude <- rownames(sobj.label)[which(Matrix::rowSums(Seurat::GetAssayData(sobj.label, slot = "data") != 0) < 2)] 
      }
      
      message(paste0("Identify ",round(ncol(sobj.label)/targetGamma)," metacells using SuperCell...\n"))
      
      SC.label <- SuperCell::SCimplify(GetAssayData(sobj.label,slot = "data"),  # normalized gene expression matrix 
                                       n.pc = n.pc,
                                       k.knn = k.knn, # number of nearest neighbors to build kNN network
                                       gamma = targetGamma, # graining level
                                       genes.use = Seurat::VariableFeatures(sobj.label)[ !Seurat::VariableFeatures(sobj.label) %in% genes_exclude])
    } else {
      message("object contain less than 5 single cells, simplification is not possible\naggregating all single-cells in one metacell.")
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
    message("Identify Metacells sequentially...\n")
  }
  
  for (label in unique(sobj[[opt$annotations]][,1])) {
    
    if(length(unique(sobj[[opt$annotations]][,1]))> 1) {
      message(paste0("Treat ",label, " cells\n"))
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
      
      sobj.label <- Seurat::FindVariableFeatures(sobj.label,nfeatures = opt$nFeatures,verbose = F) #is performed on raw counts in Seurat
      
      if(packageVersion("Seurat") >= 5) {
        genes_exclude <- rownames(sobj.label)[which(Matrix::rowSums(Seurat::GetAssayData(sobj.label, layer = "data") != 0) < 2)] 
      } else {
        genes_exclude <- rownames(sobj.label)[which(Matrix::rowSums(Seurat::GetAssayData(sobj.label, slot = "data") != 0) < 2)] 
      }
      
      
      fields <- sapply(X = colnames(sobj.label@meta.data) , 
                       FUN = function(X) {is.character(sobj.label[[X]][,1]) | is.factor(sobj.label[[X]][,1])})
      fields <- names(fields[which(fields)])
      
      message(paste0("Identify ",round(ncol(sobj.label)/targetGamma)," metacells using SuperCell...\n"))
      SC.label <- SuperCell::SCimplify(GetAssayData(sobj.label,slot = "data"),  # normalized gene expression matrix 
                                       n.pc = n.pc,
                                       k.knn = k.knn, # number of nearest neighbors to build kNN network
                                       gamma = targetGamma, # graining level
                                       genes.use = Seurat::VariableFeatures(sobj.label)[ !Seurat::VariableFeatures(sobj.label) %in% genes_exclude])
      gc(verbose = F)
    } else {
      message("object contain less than 5 single cells, simplification is not possible\naggregating all single-cells in one metacell.")
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
  
  membershipNames <- names(SC$membership)
  SC$membership <- paste0("mc",SC$membership)
  names(SC$membership) <- membershipNames
  sobj$Metacell <- SC$membership
  
  if (packageVersion("Seurat") < 5) {
    sobjMC <- AggregateExpression(sobj,
                                  return.seurat = T,
                                  group.by = "Metacell",
                                  slot = "counts", 
                                  verbose  = F)
    sobjMC@assays$RNA@data <- sobjMC@assays$RNA@counts # because AggregateExpression with slot = counts set data to log1p(aggregated counts)
  } else {
    sobjMC <- Seurat::AggregateExpression(sobj,
                                          return.seurat = T,
                                          group.by = "Metacell",
                                          verbose  = F)  #sobjMC@assays$RNA@layers$data <- sobjMC@assays$RNA@layers$counts # because AggregateExpression with slot = counts set data to log1p(aggregated counts)
    #sobjMC <- Seurat::SetAssayData(sobjMC,layer = "data",new.data = Seurat::GetAssayData(sobjMC,assay = "RNA",layer = "counts"),assay = "RNA")
    rownames(sobjMC@meta.data) <- colnames(sobjMC)
    
    
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
  
  
  
  #Store membership in sobjMC
  membershipDF <- data.frame(SC$membership)
  colnames(membershipDF) <- "membership"
  sobjMC@misc$cell_membership <- membershipDF
  
  #Store size
  sobjMC$size <- as.numeric(table(sobjMC@misc$cell_membership))
  
  
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
































