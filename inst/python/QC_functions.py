import pandas as pd
import numpy as np
import palantir
import warnings
import scipy

# adapted from https://github.com/dpeerlab/SEACells/blob/main/SEACells/evaluate.py

def get_diffusion_map(
    ad,
    low_dim_embedding='X_pca',
    n_comp = None
):
    components = pd.DataFrame(ad.obsm[low_dim_embedding]).set_index(ad.obs_names)
    
    if n_comp is None:
        n_comp = components.shape[1]
    else:
        if n_comp > components.shape[1]:
            warning(f'Number of components to use is larges than number of existing components, n_comp set to max')
            n_comp = components.shape[1]
    
    components = components.iloc[:,:n_comp]
    
    dm_res = palantir.utils.run_diffusion_maps(components, n_components=n_comp)
    dc = palantir.utils.determine_multiscale_space(dm_res, n_eigs=10)
    dc.index=components.index
  
    return(dc)

