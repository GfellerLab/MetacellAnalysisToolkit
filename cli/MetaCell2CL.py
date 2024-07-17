import numpy as np
import pandas as pd
import scanpy as sc
#import SEACells
import scanpy as sc
import os
import sys, getopt
from pathlib import Path
import metacells as mc

# from https://tanaylab.github.io/metacells-vignettes/one-pass.html

def prepare_run_metacell0_9(adata, 
                    target_metacell_size = 96, 
                    proj_name = "metacells",
                    pre_filter_cells = True,
                    EXCLUDED_GENE_NAMES = ["XIST", "MALAT1", "NEAT1"], 
                    EXCLUDED_GENE_PATTERNS = ["MT-.*"],
                    ADDITIONAL_CELLS_MASKS = None,   #name of a logical column in .obs with True for cells to discard
                    PROPERLY_SAMPLED_MIN_CELL_TOTAL = 800,
                    PROPERLY_SAMPLED_MAX_CELL_TOTAL = 20000,
                    PROPERLY_SAMPLED_MAX_EXCLUDED_GENES_FRACTION = 0.25,
                    LATERAL_GENE_NAMES = [
    "ACSM3", "ANP32B", "APOE", "AURKA", "B2M", "BIRC5", "BTG2", "CALM1", "CD63", "CD69", "CDK4",
    "CENPF", "CENPU", "CENPW", "CH17-373J23.1", "CKS1B", "CKS2", "COX4I1", "CXCR4", "DNAJB1",
    "DONSON", "DUSP1", "DUT", "EEF1A1", "EEF1B2", "EIF3E", "EMP3", "FKBP4", "FOS", "FOSB", "FTH1",
    "G0S2", "GGH", "GLTSCR2", "GMNN", "GNB2L1", "GPR183", "H2AFZ", "H3F3B", "HBM", "HIST1H1C",
    "HIST1H2AC", "HIST1H2BG", "HIST1H4C", "HLA-A", "HLA-B", "HLA-C", "HLA-DMA", "HLA-DMB",
    "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRA", "HLA-DRB1", "HLA-E", "HLA-F", "HMGA1",
    "HMGB1", "HMGB2", "HMGB3", "HMGN2", "HNRNPAB", "HSP90AA1", "HSP90AB1", "HSPA1A", "HSPA1B",
    "HSPA6", "HSPD1", "HSPE1", "HSPH1", "ID2", "IER2", "IGHA1", "IGHA2", "IGHD", "IGHG1", "IGHG2",
    "IGHG3", "IGHG4", "IGHM", "IGKC", "IGKV1-12", "IGKV1-39", "IGKV1-5", "IGKV3-15", "IGKV4-1",
    "IGLC2", "IGLC3", "IGLC6", "IGLC7", "IGLL1", "IGLL5", "IGLV2-34", "JUN", "JUNB", "KIAA0101",
    "LEPROTL1", "LGALS1", "LINC01206", "LTB", "MCM3", "MCM4", "MCM7", "MKI67", "MT2A", "MYL12A",
    "MYL6", "NASP", "NFKBIA", "NUSAP1", "PA2G4", "PCNA", "PDLIM1", "PLK3", "PPP1R15A", "PTMA",
    "PTTG1", "RAN", "RANBP1", "RGCC", "RGS1", "RGS2", "RGS3", "RP11-1143G9.4", "RP11-160E2.6",
    "RP11-53B5.1", "RP11-620J15.3", "RP5-1025A1.3", "RP5-1171I10.5", "RPS10", "RPS10-NUDT3", "RPS11",
    "RPS12", "RPS13", "RPS14", "RPS15", "RPS15A", "RPS16", "RPS17", "RPS18", "RPS19", "RPS19BP1",
    "RPS2", "RPS20", "RPS21", "RPS23", "RPS24", "RPS25", "RPS26", "RPS27", "RPS27A", "RPS27L",
    "RPS28", "RPS29", "RPS3", "RPS3A", "RPS4X", "RPS4Y1", "RPS4Y2", "RPS5", "RPS6", "RPS6KA1",
    "RPS6KA2", "RPS6KA2-AS1", "RPS6KA3", "RPS6KA4", "RPS6KA5", "RPS6KA6", "RPS6KB1", "RPS6KB2",
    "RPS6KC1", "RPS6KL1", "RPS7", "RPS8", "RPS9", "RPSA", "RRM2", "SMC4", "SRGN", "SRSF7", "STMN1",
    "TK1", "TMSB4X", "TOP2A", "TPX2", "TSC22D3", "TUBA1A", "TUBA1B", "TUBB", "TUBB4B", "TXN", "TYMS",
    "UBA52", "UBC", "UBE2C", "UHRF1", "YBX1", "YPEL5", "ZFP36", "ZWINT"
],
                    LATERAL_GENE_PATTERNS = ["RP[LS].*"],
                    NOISY_GENE_NAMES = [
    "CCL3", "CCL4", "CCL5", "CXCL8", "DUSP1", "FOS", "G0S2", "HBB", "HIST1H4C", "IER2", "IGKC",
    "IGLC2", "JUN", "JUNB", "KLRB1", "MT2A", "RPS26", "RPS4Y1", "TRBC1", "TUBA1B", "TUBB"
],
                    seed = 123456
):
                  
  mc.ut.set_name(adata, proj_name)
  # exclude genes based on prior knowledge + automatically
  adata.X.sort_indices()
  
  mc.pl.exclude_genes(
    adata,
    excluded_gene_names=EXCLUDED_GENE_NAMES, 
    excluded_gene_patterns=EXCLUDED_GENE_PATTERNS,
    random_seed=seed
  )
  
  if pre_filter_cells:
      print("Filter cells using standard MetaCell pipeline...")
      # compute number of UMIs in excluded genes 
      mc.tl.compute_excluded_gene_umis(adata)
      
      if ADDITIONAL_CELLS_MASKS is not None:
          ADDITIONAL_CELLS_MASKS = adata.obs[ADDITIONAL_CELLS_MASKS]
  
      # exclude cells based on totals UMIs and number of UMIs in excluded genes + additional cell annotation provided by the user (additional_cells_masks)
      mc.pl.exclude_cells(
      adata,
      properly_sampled_min_cell_total=PROPERLY_SAMPLED_MIN_CELL_TOTAL,
      properly_sampled_max_cell_total=PROPERLY_SAMPLED_MAX_CELL_TOTAL,
      properly_sampled_max_excluded_genes_fraction=PROPERLY_SAMPLED_MAX_EXCLUDED_GENES_FRACTION,
      additional_cells_masks=ADDITIONAL_CELLS_MASKS
      )
  else:
      mc.pl.exclude_cells(
      adata,
      properly_sampled_min_cell_total=None,
      properly_sampled_max_cell_total=None,
      properly_sampled_max_excluded_genes_fraction=None,
      additional_cells_masks=None
      )
      
      
  
  clean = mc.pl.extract_clean_data(adata, name="hca_bm.one-pass.clean")
  mc.ut.top_level(clean) 
  print(f"Clean: {clean.n_obs} cells, {clean.n_vars} genes")
  
  # Define lateral and noisy genes 
  mc.pl.mark_lateral_genes(
      clean,
      lateral_gene_names=LATERAL_GENE_NAMES,
      lateral_gene_patterns=LATERAL_GENE_PATTERNS,
  )
  mc.pl.mark_noisy_genes(clean, noisy_gene_names=NOISY_GENE_NAMES)
  
  return(clean)



def main(argv):
    # Default setting
    outdir = '.'
    pre_filter_cells = False
    gamma = 75
    annotations = None
    min_metacells = 1
    output = "adata"
    yml_config = None
    #default setting for MetaCell2
    ADDITIONAL_CELLS_MASKS = None
    seed = 123456
    
    try:
        opts, args = getopt.getopt(argv,"i:o:g:y:ps:a:m:",["input_file=","outdir=","gamma=","yml_config=","pre_filter_cells=","output=","annotations=","min_metacells="])
    except getopt.GetoptError:
        print('MetaCell2CL.py -i <input_file> -o <outdir> -g <gamma> -y <yml_config> -p <pre_filter_cells> -s <output> -a <annotations> -m <min_metacells>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('MetaCell2CL.py -i <input_file> -o <outdir> -g <gamma> -y <yml_config> -p <pre_filter_cells> -s <output> -a <annotations> -m <min_metacells>')
            sys.exit()
        elif opt in ("-i", "--input_file"):
            input_file = arg
        elif opt in ("-o", "--outdir"):
            outdir = arg
        elif opt in ("-g", "--gamma"):
            gamma = float(arg)
        elif opt in ("-y", "--yml_config"):
            yml_config = arg
        elif opt in ("-p", "--pre_filter_cells"):
            pre_filter_cells = True
        elif opt in ("-s", "--output"):
            output = arg
        elif opt in ("-a", "--annotations"):
            annotations = arg
        elif opt in ("-m", "--min_metacells"):
            min_metacells = int(arg)

    print('input is "', input_file)
    print('Output dir is "', outdir)
    print('gamma is "', gamma)
    print('Pre filter cells is',pre_filter_cells)

    
    os.makedirs(outdir,exist_ok = True)
    
    if input_file.endswith(".rds"):
        try:
            import logging
            import anndata2ri
            import rpy2.rinterface_lib.callbacks
            import rpy2.rinterface_lib.embedded
            import rpy2.robjects as ro
        
            rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)
        except ModuleNotFoundError as e:
            raise OptionalDependencyNotInstalled(e)
    
        try:
            ro.r("library(Seurat)")
            ro.r("library(scater)")
        except rpy2.rinterface_lib.embedded.RRuntimeError as ex:
            RLibraryNotFound(ex)
        
        anndata2ri.activate()
    
        ro.r(f'sobj <- readRDS("{input_file}")')
        adata = ro.r("as.SingleCellExperiment(sobj)")
        if "logcounts" in adata.layers:
            del adata.layers['logcounts'] # we only load raw counts stored in adata.X when using as.SingleCellExperiment and annadata2ri
    
        anndata2ri.deactivate()
    
    else:
        adata = sc.read_h5ad(input_file)
        
        if adata.raw is not None:
            adata = adata.raw.to_adata()
            # adata.X = adata.raw.X # we only load raw counts, We always normalize .X prior to compute PCA if prePro is asked or reduction_key absent  
            # del adata.raw

    # The dtype of X is no longer set to float32 in scampy. 
    # While anndata2ri produces float64, the majority of h5ad objects available online are float32.
    # We choose to set the type to float32
    
    adata.X = adata.X.astype("float32")
    
    if yml_config is not None:
        print("Using MetaCell2 with parameters from"+ yml_config)
        with open(yml_config) as file:
            mc2_params = yaml.safe_load(file)
            
        if mc2_params["ADDITIONAL_CELLS_MASKS"] == "None":
            mc2_params["ADDITIONAL_CELLS_MASKS"] = None
        
        clean = prepare_run_metacell0_9(
          adata,
          target_metacell_size = gamma,
          pre_filter_cells = pre_filter_cells,
          EXCLUDED_GENE_NAMES = mc2_params["EXCLUDED_GENE_NAMES"], 
          EXCLUDED_GENE_PATTERNS = mc2_params["EXCLUDED_GENE_PATTERNS"],
          ADDITIONAL_CELLS_MASKS = mc2_params["ADDITIONAL_CELLS_MASKS"],
          PROPERLY_SAMPLED_MIN_CELL_TOTAL =  mc2_params["PROPERLY_SAMPLED_MIN_CELL_TOTAL"],
          PROPERLY_SAMPLED_MAX_CELL_TOTAL = mc2_params["PROPERLY_SAMPLED_MAX_CELL_TOTAL"],
          PROPERLY_SAMPLED_MAX_EXCLUDED_GENES_FRACTION = mc2_params["PROPERLY_SAMPLED_MAX_EXCLUDED_GENES_FRACTION"],
          LATERAL_GENE_NAMES = mc2_params["LATERAL_GENE_NAMES"],
          LATERAL_GENE_PATTERNS = mc2_params["LATERAL_GENE_PATTERNS"],
          NOISY_GENE_NAMES = mc2_params["NOISY_GENE_NAMES"],
          seed = mc2_params["seed"])
          
    else:
        clean = prepare_run_metacell0_9(
          adata,
          target_metacell_size = gamma,
          pre_filter_cells = pre_filter_cells)
    
  
    print("Identify"+str(clean.n_obs/gamma)+"metacells using MetaCell2")  
    # Adapt number of parallel piles to the available memory 
    max_parallel_piles = mc.pl.guess_max_parallel_piles(clean)
    mc.pl.set_max_parallel_piles(max_parallel_piles)
    mc.pl.divide_and_conquer_pipeline(
      clean,
      target_metacell_size = gamma,
      random_seed = seed)
  
    # Aggregate metacells 
    adata_mc = mc.pl.collect_metacells(clean, name='metacells', random_seed = seed)
    #Store membership
    clean.obs['membership'] = [i+1 if i >= 0 else np.nan for i in clean.obs["metacell"]]
    # Store it in adata_mc
    adata_mc.uns["cell_membership"] = clean.obs[["membership","metacell_name"]]
    # Rename metacells size column
    adata_mc.obs = adata_mc.obs.rename(columns={'grouped':'size'})
    
    print("Assign metadata to metacells and compute purities...")
    #We use SEACells function for this
    from SEACells.evaluate import compute_celltype_purity
    clean.obs["SEACell"] = clean.obs['metacell_name']
    for c in clean.obs.select_dtypes(exclude='number').columns:
        if c not in [ 'properly_sampled_cell', 'excluded_cell','dissolved', 'rare_cell', 'metacell_name',"SEACells"]:
            metacell_purity = compute_celltype_purity(clean, c)
            adata_mc.obs = adata_mc.obs.join(metacell_purity)



    if output == "adata":
        
        print("Save results as adata...")
        adata_out = outdir+"/mc_adata.h5ad"
        adata_mc.write_h5ad(adata_out)


    if output == "seurat":
        
        print("Save results as seurat...")
        seurat_out = outdir+"/mc_Seurat.rds"
    
        try:
            import logging
            from scipy import sparse 
            import anndata2ri
            import rpy2.rinterface_lib.callbacks
            import rpy2.rinterface_lib.embedded
            import rpy2.robjects as ro
        
            rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)
    
        except ModuleNotFoundError as e:
            raise OptionalDependencyNotInstalled(e)
      
        try:
            ro.r("library(Seurat)")
            ro.r("library(scater)")
        except rpy2.rinterface_lib.embedded.RRuntimeError as ex:
            RLibraryNotFound(ex)
        
        anndata2ri.activate()
    
        if sparse.issparse(adata_mc.X):
            if not adata_mc.X.has_sorted_indices:
                adata_mc.X.sort_indices()
        
        #Would be needed if other layers are summarized    
        #for key in adataMC.layers:
        #    if sparse.issparse(adataMC.layers[key]):
        #        if not adataMC.layers[key].has_sorted_indices:
        #            adataMC.layers[key].sort_indices()
                
        ro.globalenv["adata_mc"] = adata_mc
        ro.globalenv["membership"] = adata_mc.uns['cell_membership']
    
        ro.r('''
        sobj.mc <- CreateSeuratObject(counts = assay(adata_mc),meta.data = data.frame(colData(adata_mc)))
        sobj.mc@misc$cell_membership <- membership  
        ''')
    
        ro.r(f'saveRDS(sobj.mc, file="{seurat_out}")')
    
        anndata2ri.deactivate()



if __name__ == "__main__":
    main(sys.argv[1:])
    






