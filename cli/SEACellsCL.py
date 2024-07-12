import numpy as np
import pandas as pd
import scanpy as sc
import SEACells
import scanpy as sc
import os
import sys, getopt
from pathlib import Path

def main(argv):
    # Default setting
    outdir = '.'
    dim_str = "1:50"
    gamma = 75
    n_features = 2000
    reduction_key = None
    annotations = None
    min_metacells = 1
    output = "adata"
    #default setting for SEACells
    n_waypoint_eigs = 10 # Number of eigenvalues to consider when initializing metacells
    min_iter=10
    max_iter=100
    k_knn = 15
    
    try:
        opts, args = getopt.getopt(argv,"i:o:r:n:f:k:g:s:a:m:",["input_file=","outdir=","reduction_key=","dims=","n_features=","k_knn","gamma=","output=","annotations=","min_metacells="])
    except getopt.GetoptError:
        print('SEACellsCL.py -i <input_file> -o <outdir> -r <reduction_key> -n <dims> -f <n_features> -k <k_knn> -g <gamma> -s <output> -a <annotations> -m <min_metacells>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('SEACellsCL.py -i <input_file> -o <outdir> -r <reduction_key> -n <dims> -f <n_features> -k <k_knn> -g <gamma> -s <output> -a <annotations> -m <min_metacells>')
            sys.exit()
        elif opt in ("-i", "--input_file"):
            input_file = arg
        elif opt in ("-o", "--outdir"):
            outdir = arg
        elif opt in ("-r", "--reduction_key"):
            reduction_key= arg
            print('reduction_key is"', reduction_key)
        elif opt in ("-n", "--dims"):
            dim_str= arg
        elif opt in ("-f", "--n_features"):
            n_features= int(arg)
        elif opt in ("-k", "--k_knn"):
            k_knn= int(arg)
        elif opt in ("-g", "--gamma"):
            gamma = float(arg)
        elif opt in ("-s", "--output"):
            output = arg
        elif opt in ("-a", "--annotations"):
            annotations = arg
        elif opt in ("-m", "--min_metacells"):
            min_metacells = int(arg)

    print('input is "', input_file)
    print('Output dir is "', outdir)
    print('gamma is "', gamma)
    print('dims are "', dim_str)
    print('k knn is"', k_knn)
    
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
            del adata.layers['logcounts']  # we only load raw counts stored in adata.X when using as.SingleCellExperiment and annadata2ri
    
        anndata2ri.deactivate()
    
    else:
        adata = sc.read_h5ad(input_file)
        
        if adata.raw is not None:
            adata.X = adata.raw.X  # we only load raw counts, We always normalize .X prior to compute PCA if reduction_key absent or not in the object  
            del adata.raw

    # The dtype of X is no longer set to float32 in scampy. 
    # While anndata2ri produces float64, the majority of h5ad objects available online are float32.
    # We choose to set the type to float32
    
    adata.X = adata.X.astype("float32")
    
    dim_str_list = dim_str.split(":") # range input is interesting when using SEACells for ATAC data for which 1 component is often discarded
    # 1 to given components are used when only one number is given to dim
    if (len(dim_str_list)<2):
        dim_str_list += "1"
        dim_str_list.reverse()
    
    # Copy the counts to ".raw" attribute of the anndata since it is necessary for downstream analysis
    # This step should be performed after filtering
    raw_ad = sc.AnnData(adata.X)
    raw_ad.obs_names, raw_ad.var_names = adata.obs_names, adata.var_names
    adata.raw = raw_ad

    if annotations is None:
        annotations = "SEACell_batch"
        adata.obs["SEACell_batch"] = "all_cells"

    for anno in adata.obs[annotations].unique():
        adata_label = adata[adata.obs[annotations] == anno,]
        n_SEACells = round(len(adata_label)/gamma)
        
        if n_SEACells < min_metacells:
            n_SEACells = min_metacells
            
        if n_SEACells == 1:
            adata_label.obs['SEACell'] = "SEACell-1_"+ anno
            if anno == adata.obs[annotations].unique()[0]:
                seacells_res = adata_label.obs["SEACell"]
            else:
                seacells_res = pd.concat([seacells_res,adata_label.obs["SEACell"]])
            continue
        
        print("Identify "+ str(n_SEACells) + " metacells using SEACells...")
        
        if (reduction_key is None or reduction_key not in adata_label.obsm.keys()):
            print("Preprocess the data...")
            print("Normalize cells and compute highly variable genes...")
            sc.pp.normalize_per_cell(adata_label)
            sc.pp.log1p(adata_label)
            sc.pp.highly_variable_genes(adata_label, n_top_genes=n_features)
        
            print("Compute principal components")
            sc.tl.pca(adata_label, n_comps=int(dim_str_list[1]), use_highly_variable=True)
            reduction_key = "X_pca"
        else:
            print("using pre-computed dimension reduction "+ reduction_key)
        
        build_kernel_on = reduction_key
        
        if int(dim_str_list[1]) > len(adata_label.obsm[reduction_key][0]):
            print("number of PCs requested superior to PCs available setting last PC number to the number of available PCs")
            dim_str_list[1] = len(adata_label.obsm[reduction_key][0])

        adata_label.obsm[build_kernel_on] = adata_label.obsm[build_kernel_on][:,range(int(dim_str_list[0])-1, int(dim_str_list[1]))]
        
    
        min_metacells = min(min_metacells,adata_label.n_obs)
    
        
        
        if n_SEACells < n_waypoint_eigs:
            n_waypoint_eigs_label = n_SEACells
        else:
            n_waypoint_eigs_label = n_waypoint_eigs
      
    
        model = SEACells.core.SEACells(adata_label, 
        build_kernel_on=build_kernel_on, 
        n_SEACells=n_SEACells, 
        n_neighbors = k_knn,
        n_waypoint_eigs=n_waypoint_eigs_label,
        convergence_epsilon = 1e-5)
        
        model.construct_kernel_matrix()
        # M = model.kernel_matrix
        model.initialize_archetypes()
        #SEACells.plot.plot_initialization(ad, model)
        model.fit(min_iter=min_iter, max_iter=max_iter)
        #model.plot_convergence()
        adata_label.obs['SEACell'] = adata_label.obs['SEACell'] +"_"+ anno 
        
        if anno == adata.obs[annotations].unique()[0]:
            seacells_res = adata_label.obs["SEACell"]
        else:
            seacells_res = pd.concat([seacells_res,adata_label.obs["SEACell"]])
        

    adata.obs['SEACell'] = seacells_res.reindex(adata.obs_names)
    adata_mc = SEACells.core.summarize_by_SEACell(adata, SEACells_label='SEACell', summarize_layer='raw')
        
    # Store metacells size
    label_df = adata.obs[['SEACell']].reset_index()
    adata_mc.obs = adata_mc.obs.join(pd.DataFrame(label_df.groupby('SEACell').count().iloc[:, 0]).rename(columns={'index':'size'}))
        
    #Store membership
    d = {x: i+1 for i, x in enumerate(adata_mc.obs_names)}
    # make a membership (as in SCimplify() from SuperCell) vector
    adata.obs['membership'] = [d[x] for x in seacells_res]
    # Store it in adata_mc
    adata_mc.uns['cell_membership'] = adata.obs[["membership","SEACell"]]


    print("Assign metadata to metacells and compute purities...")
    for c in adata.obs.select_dtypes(exclude='number').columns:
        if c != "SEACell" and c != "SEACell_batch":
            SEACell_purity = SEACells.evaluate.compute_celltype_purity(adata, c)
            adata_mc.obs = adata_mc.obs.join(SEACell_purity)

    # writing membership dataframe
    # adata.obs[["SEACell","membership"]].to_csv(outdir + "SEACells_memberships.csv",index_label= "cell")
    
    #iterative merging of adata per cell annotation
    #loop end

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
    






