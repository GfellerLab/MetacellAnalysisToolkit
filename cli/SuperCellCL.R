# library(Seurat)
library(getopt)
#library(doParallel)
# library(SuperCell)

spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'input',  'i', 1, "character", "REQUIRED : Path to seurat object (.rds) or anndata (.h5ad) object",
  'outdir',     'o', 1, "character", 'Outdir path (default ./)',
  'isNorm', 'd', 0, "logical", "Whether the data are already normalised (default FALSE and log normalize data using Seurat)",
  "n.pc", "n", 1, "numeric", "Number of principal components to use for construction of single-cell kNN network (default 50)",
  "nFeatures", "f", 1, "numeric", "number of highly variable features to use (default 2000), selected with Seurat",
  "k.knn", "k", 1, "numeric", "k for the knn used in the wnn used for metacell identification (default 30)",
  "gamma", "g", 1, "numeric", "gamma used for metacell identification (default 10)",
  "annotations", "a", 1, "character", "Identify metacell per label regarding a gigven annotation (column name in the metadata of the object)",
  "cores", "l", 1, "numeric", "number of cores to use for parallel identifications of metacells in each label",
  "minMetacells", "m", 1, "numeric", "minimum number of metacells to identify (usefull and annotations is given)",
  "output", "s", 1, "character", "how to save the SuperCell results ('seurat' : Seurat object (.rds), 'adata' : anndata object (.h5ad), 'SC' : SuperCell object (SC.rds))"
), byrow=TRUE, ncol=5)

opt = getopt(spec)

if ( !is.null(opt$help) | (is.null(opt$input) & is.null(opt$output))) {
  cat("Command line implemention of SuperCell workflow:")
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

#if help was asked, print a friendly message
#and exit with a non-zero error code
#test
# setwd("~/Documents/reviewTutorial/")
# opt <- list()
# opt$input <- "data/HCLA/datasets/Banovich_Kropski_2020/sc_adata.h5ad"
# #opt$doNorm <- T
# opt$n.pc <- 30
# opt$k.knn <- 30
# opt$minMetacells <- 50
# opt$outdir = "data/HCLA/datasets/Banovich_Kropski_2020"
# opt$gamma = 50
# opt$nFeature = 2000
# opt$output <- "adata"
# opt$annotations <- "sample"
# opt$cores <- 6
# print(opt)




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
    sobj <- Seurat::CreateSeuratObject(counts = countMatrix,meta.data = adata$obs)
    remove(countMatrix)
  } else { # we assume norm counts are in .X and raw counts in .raw.X
    normMatrix <-  Matrix::t(adata$X)
    colnames(normMatrix) <- adata$obs_names
    rownames(normMatrix) <- adata$var_names
    sobj <- Seurat::CreateSeuratObject(counts = normMatrix,meta.data = adata$obs)
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
  print("Normalize data...")
  sobj <- Seurat::NormalizeData(sobj,verbose = F)
}

fields <- sapply(X = colnames(sobj@meta.data) , 
                 FUN = function(X) {is.character(sobj[[X]][,1]) | is.factor(sobj[[X]][,1])})
fields <- names(fields[which(fields)])

print("Identify Metacells...")


if (opt$cores > 1 & !is.null(opt$annotations)) {
  
  SCs <- foreach::foreach(sobj.label = Seurat::SplitObject(sobj,split.by = opt$annotations)) %dopar% {
    
    #print(paste0("Treat ",label, " cells"))
    #sobj.label <- sobj[,sobj[[opt$annotations]][,1] == label]
    
    # adapt parameters regarding number of single cells (occurs mainly when an annotation is given)
    if (is.null(opt$minMetacells)) {
      opt$minMetacells <- 1
    }
    minMetacells <- min(ncol(sobj.label),opt$minMetacells)
    k.knn <- min(opt$k.knn, ncol(sobj.label)-1)
    n.pc <- min(opt$n.pc, ncol(sobj.label)-1)
    
    
    
    targetGamma <- min(ncol(sobj.label)/minMetacells,opt$gamma)
    

    
    
    sobj.label <- Seurat::FindVariableFeatures(sobj.label,nfeatures = opt$nFeatures,verbose = F) #is performed on raw counts (as in Seurat workflow) only if is.norm = F 
    
    
    
    print(paste0("Identify ",round(ncol(sobj.label)/targetGamma)," metacells using SuperCell..."))
    
    SC.label <- SuperCell::SCimplify(Seurat::GetAssayData(sobj.label,slot = "data"),  # normalized gene expression matrix 
                                     n.pc = n.pc,
                                     k.knn = k.knn, # number of nearest neighbors to build kNN network
                                     gamma = targetGamma, # graining level
                                     genes.use = Seurat::VariableFeatures(sobj.label))
    
    #SCs[[label]] <- SC.label
  }
  
  parallel::stopCluster(cluster)
  gc()
  
} else {
  print("Identify Metacells sequentially...")
  
  SCs <- list()
  if (is.null(opt$annotations)) {
    opt$annotations <- "orig.ident"
  }
  
  for (label in unique(sobj[[opt$annotations]][,1])) {
    
    print(paste0("Treat ",label, " cells"))
    sobj.label <- sobj[,sobj[[opt$annotations]][,1] == label]
    
    # adapt parameters regarding number of single cells (occurs mainly when an annotation is given)
    if (is.null(opt$minMetacells)) {
      opt$minMetacells <- 1
    }
    minMetacells <- min(ncol(sobj.label),opt$minMetacells)
    k.knn <- min(opt$k.knn, ncol(sobj.label)-1)
    n.pc <- min(opt$n.pc, ncol(sobj.label)-1)
    
    
    
    targetGamma <- min(ncol(sobj.label)/minMetacells,opt$gamma)
    
    sobj.label <- Seurat::FindVariableFeatures(sobj.label,nfeatures = opt$nFeatures,verbose = F) #is performed on raw counts in Seurat
    
    
    fields <- sapply(X = colnames(sobj.label@meta.data) , 
                     FUN = function(X) {is.character(sobj.label[[X]][,1]) | is.factor(sobj.label[[X]][,1])})
    fields <- names(fields[which(fields)])
    
    print(paste0("Identify ",round(ncol(sobj.label)/targetGamma)," metacells using SuperCell..."))
    if (!exists("normMatrix")) {
      SC.label <- SuperCell::SCimplify(Seurat::GetAssayData(sobj.label,slot = "data"),  # normalized gene expression matrix 
                                       n.pc = n.pc,
                                       k.knn = k.knn, # number of nearest neighbors to build kNN network
                                       gamma = targetGamma, # graining level
                                       genes.use = Seurat::VariableFeatures(sobj.label))
    } else {
      SC.label <- SuperCell::SCimplify(normMatrix,  # normalized gene expression matrix in case of normalized anndata object as input
                                       n.pc = n.pc,
                                       k.knn = k.knn, # number of nearest neighbors to build kNN network
                                       gamma = targetGamma, # graining level
                                       genes.use = Seurat::VariableFeatures(sobj.label))
      remove(normMatrix)
      gc(verbose = F)
    }
    
    SCs[[label]] <- SC.label
  }
}

SC <- SuperCell::supercell_merge(SCs)

if (opt$output != "SC") {
  
  if (opt$isNorm) {
    sobj <- Seurat::CreateSeuratObject(counts = countMatrix,meta.data = sobj@meta.data)
    gc()
  }
  
  sobj$Metacell <- SC$membership
  sobjMC <- Seurat::AggregateExpression(sobj,
                                        return.seurat = T,
                                        group.by = "Metacell",
                                        slot = "counts", verbose  = F)
  sobjMC@assays$RNA@data <- sobjMC@assays$RNA@counts # because AggregateExpression with slot = counts set data to log1p(aggregated counts)
  gc()
  
  
  print("Assign metadata to metacells and compute purities...")
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
    adataMC <- anndata::AnnData(X = Matrix::t(Seurat::GetAssayData(object = sobjMC,slot = "counts")),
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
































