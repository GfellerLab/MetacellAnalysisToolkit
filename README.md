# MetacellAnalysisToolkit (MATK)

Toolkit for metacell analysis. It consists of the `MATK` command line tool to easily identify metacells with either [SEACells](https://github.com/dpeerlab/SEACells), [SuperCell](https://github.com/GfellerLab/SuperCell) or [MetaCell2](https://github.com/tanaylab/metacells/tree/master) with a joined R package for metacell quality control and visualization.

## 1. Installation

### 1.1 Clone the GitHub repository and move to the MetacellAnalysisToolkit directory

### 1.2 Create the conda environment

Then you need to create the conda environment that contains most of python and R packages useful for metacell analyses. You need to have a conda installer such as [Miniconda](https://docs.conda.io/projects/miniconda/en/latest/).

    conda env create -f env/MetacellAnalysisToolkit_env.yml

Alternatively you can use [mamba](https://github.com/conda-forge/miniforge) for a faster and lighter installation.

    mamba env create -f env/MetacellAnalysisToolkit_env.yml

### 1.3 Install additional R packages

Then you have to install in this environment additional required R packages not available through conda.

    conda activate MetacellAnalysisToolkit
    Rscript env/install.R

#### Seurat v5 compatibility

This toolkit has been developped under seurat version 4 which is the Seurat version installed with the MetacellAnalysisToolkit environment.
We recommand using Seurat v4 and this environment to use MATK. However, this toolkit is also compatible with Seurat v5 installed by using the v4 assay [option](https://satijalab.org/seurat/articles/seurat5_essential_commands#create-seurat-or-assay-objects).

### 1.4 Make command line scripts executable

    chmod a+x cli/MATK
    chmod a+x cli/SuperCellCL.R 
    chmod a+x cli/SEACellsCL.py
    chmod a+x cli/MetaCell2CL.py

### 1.5 Configure PATH

If you want, you can finally add the value of the path to the `cli` directory (of this repository) to your PATH environment variable so that you can use the MATK command line tool directly. On Linux, using bash, You can do this by adding this line to your `~/.bashrc` (or `~/.bash_profile` on macOS):

    export PATH="/path/to/MetacellAnalysisToolkit/cli/:$PATH"

Don't forget to source your `~/.bashrc` (or `~/.bash_profile` on macOS) after.

### Use of MATK within a Docker container

We also provide a Docker file to build an environment with all the requirements to run MATK. You can build the docker environment using the following command line:

    docker build -t matk:v1.1 -f env/Dockerfile_MATK .
    
On MAC, if you are encountering issues, try the following command line:

    docker build --platform linux/amd64 -t matk:v1.1 -f env/Dockerfile_MATK .

You can also pull our prebuilt image using: 

    docker pull agabriel/matk:v1.1

To run MATK on a test dataset (downloaded in section 2) within this docker container with docker or singularity please refer to section 3.5.

Note that the container corresponding to the dockerfile `env/Dockerfile_MATK` is based on Seurat V5, if you want to use Seurat V4, use `env/Dockerfile_MATK_SeuratV4` or use the following prebuilt image: `agabriel/matk:v1.0`

## 2. Download test data

MATK takes as input/output either an Anndata .h5ad objects or Seurat .rds object.

### 2.1 CD34+ scRNA-seq dataset (6,900 cells) from Dana's Peer lab (.h5ad file).

    wget https://zenodo.org/records/6383269/files/cd34_multiome_rna.h5ad?download=1 -O data/cd34_multiome_rna.h5ad

### 2.2 PBMC scRNA-seq (6,900 cells) dataset from scanpy datasets

Here we use short python and R scripts to get a .h5ad and a .rds object

    python get_data/get_PBMC_dataset.py
    Rscript get_data/get_PBMC_rds.R

## 3. Usage

Using MATK tool you can easily identify metacells with either SEACells, SuperCell or MetaCell2 using various common and method-specific options

### 3.1 Print help

    $MATK -h
    usage: /path/to/MetacellAnalysisToolkit/cli/MATK options

    Constructing metacell from single cell data with SEACells (0.3.3) 'MetaCell2 (0.9.0) or SuperCell (1.0)
    Expect a filtered (low quality cells removed) Seurat or Anndata object  

    1 - Identifying metacells, 
    2 - aggregating counts data per metacell (summing raw counts)
    3 - assigning metadata to metacells and computing purities (Assigning metacells to the most aboundant label)

    OPTIONS:
       -h     Show this message
   
       -t     tool, either 'SEACells', 'MetaCell' or 'SuperCell' 
   
       -i     input_file, either an Anndata object file '.h5ad' or a Seurat object file '.rds' file
   
       -o     outdir, output directory (default ./)

       -n     dims, number of principal components to use (only for SEACells and SuperCell, default 50) 
   
       -f     n_features, number of highly variable genes use to compute the initial PCA (only for SEACells and SuperCell, default 2000) 
   
       -k     k_knn, number of neighbors to construct the knn graph (only for SEACells and SuperCell, default 30)
   
       -g     gamma,        graining level of data 
              Proportion of number of single cells in the initial dataset to the number of metacells in the final dataset
              When using MetaCell this correspond to a target gamma (obtained gamma slightly lower)
          
       -s     output, desired metacell file format in output, either 'adata' for a h5ad file or a 'seurat' for a rds file. 
              Output file name will be  'mc_'{output_format}. 
          
       -d     normalized data in input (only for SuperCell). ADD -d to specify that data are already normalized in the data slot of the Seurat object or in .X for a adata object (default FALSE).
              Note that in this case raw count data have to be provided in the count slot for a Seurat object or in .raw.X for an anndata object. 
   
       -r     reduction_key (only for SEACells, default none and a PCA reduction is computed using standard scanpy workflow and stored in "X_pca")
   
       -y     yaml_file (only for MetaCell2, default None and use default options and gene lists)
   
       -a     annotation, to make supervised metacells according to an annotation present in the metadata of the object (only for SEACells and SuperCell, default none)
   
       -l     cores, number of cores to use for parallel processing if an annotation is profided (only for SuperCell)


### 3.2 Metacell identification on Cd34+ cells using SuperCell

Here we identify metacells from the h5ad file of CD34+ cells using SuperCell and save the results in a h5ad file. We use 50 principal components, 30 neighbors for the knn and a graining level of 75.

    MATK -t SuperCell -i  data/cd34_multiome_rna.h5ad -o MATK_output/SuperCell/cd34/ -n 50 -f 2000 -k 30 -g 75 -s adata

### 3.2 Metacell identification on PBMCs using SEACells

Here we identify metacells from the rds file of PBMCs using SEACells and save the results in a rds file. We use here a graining level of 50.

    MATK -t SEACells -i data/pbmc.rds -o MATK_output/SEACells/pbmc/ -n 50 -f 2000 -k 30 -g 50 -s seurat

### 3.3 Metacell identification on PBMCs using MetaCell

Here we identify metacells from the rds file of PBMCs using MetaCell (v2 python version) and save the results in a h5ad file.

    MATK -t MetaCell -i data/pbmc.rds -o MATK_output/MetaCell/pbmc/ -g 50 -s seurat

MetaCell does not use a knn graph from PCA based the highly variable genes but has its own parameters (including different gene lists) you can set using a yaml config file and the `-y` argument. You have an example of a such yaml file [here]((/cli/config/MetaCell2_config.yml)) containing the default settings of MetaCell proposed by the authors.

### 3.4 Supervised Metacell identification

You can identify metacells according to a given annotation (e.g. cell types, samples) present in the metadata of the object using the `-a` argument. You can specify a minimum number of metacells to identify (per annotation) using the `-m` argument.

    #Using SEACells
    python cli/SEACellsCL.py -i data/cd34_multiome_rna.h5ad -o testCLI/SEACells_per_celltype_min_5_MC/cd34_multiome_rna/input_raw_adata/ -a celltype -m 5 -n 50 -f 2000 -k 30 -g 75 -s adata

With SuperCell it is possible to use parallel processing using the `-l` argument which gives the number of cores to use.

    #SuperCell parallel metacell identification in each cell type
    Rscript cli/SuperCellCL.R -i data/cd34_multiome_rna.adata -o testCLI/SuperCell_per_celltype/cd34_multiome_rna/input_raw_adata/ -n 50 -f 2000 -k 30 -g 75 -s adata -a celltype -l 6
    
### 3.5 Run MATK within the docker container.

To run MATK on the CD34 dataset within the docker container, use the following command line:

    docker run --rm -v $(pwd):/workspace -v $(pwd):/workspace agabriel/matk:v1.1 MATK -t SuperCell -i /workspace/data/cd34_multiome_rna.h5ad -o /workspace/MATK_output/SuperCell/cd34/ -n 50 -f 2000 -k 30 -g 75 -s adata

You can also use the container with singularity, for example to use MATK on a cluster :

    singularity pull docker://agabriel/matk:v1.1 
    singularity run --bind $(pwd) matk_v1.1.sif MATK -t SuperCell -i  data/cd34_multiome_rna.h5ad -o MATK_output/SuperCell/cd34/ -n 50 -f 2000 -k 30 -g 75 -s adata
    
## Quality control visualization

[Perform quality controls on metacells using MetacellAnalysisToolkit R package](https://github.com/GfellerLab/MetacellAnalysisToolkit/blob/main/vignettes/MetacellAnalysisToolkit_vignette.md)

## Advanced analysis

-   [Analysis](/examples/HLCA_core_atlas.md) of the core HLCA atlas comprising 500'000 cells at the metacell level using SuperCell in command line and Seurat-rpca integration.
-   [Supervised Analysis](/examples/HLCA_core_atlas_supervised.md) of the core HLCA atlas comprising 500'000 cells at the metacell level using SuperCell in command line and [STACAS](https://github.com/carmonalab/STACAS) integration.

    

