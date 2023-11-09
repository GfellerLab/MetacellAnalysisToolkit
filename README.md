# MetacellAnalysisToolkit (MCAT)

Toolkit for metacell analysis. It consists of the `MCAT` command line tool to easily identify metacells with either SEACells, SuperCell or MetaCell2 with a joined R package for metacell quality control and visualization.

## 1. Installation

### 1.1 Clone the GitHub repository and move to the MetacellAnalysisToolkit directory

### 1.2 Create the conda environment

Then you need to create the conda environment that contains most of python and R packages useful for metacell analyses

    conda env create -f MetacellToolkit_env.yml

### 1.3 Install additional R packages

Then you have to install in this environment additional required R packages not available through conda

    conda activate MetacellToolkit
    Rscript install.R

### 1. Make command line script executable

    chmod a+x cli/MCAT
    chmod a+x cli/SuperCellCL.R 
    chmod a+x cli/SEACells.py
    chmod a+x cli/MetaCell2CL.py

### 1.5 Configure PATH

If you want you can finally add the value of the path to the `cli` directory of this repository to your PATH environment variable so that you can use the MCAT command line tool directly. On Linux, using bash, You can do this by adding this line to your \~/.bashrc:

    export PATH="/home/your_account/path_to_cli_dir/:$PATH"

## 2. Download test data

MCAT takes as input/output either an Anndata .h5ad objects or Seurat .rds object.

### 2.1 CD34+ scRNA-seq dataset (6,900 cells) from Dana's Peer lab (.h5ad file).

    mkdir data/
    wget https://dp-lab-data-public.s3.amazonaws.com/SEACells-multiome/cd34_multiome_rna.h5ad -O data/cd34_multiome_rna.h5ad

### 2.2 PBMC scRNA-seq (6,900 cells) dataset from scanpy datasets

Here we use short python and R scripts to get a .h5ad and a .rds object

    python data/get_PBMC_dataset.py
    Rscript data/get_PBMC_rds.R

## 3. Usage

Using MCAT tool you can easily identify metacells with either SEACells, SuperCell or MetaCell2 using various common and method-specific options

### 3.1 Print help

    $MetacellToolkit.sh -h
    usage: /home/leonard/Documents/reviewTutorial/MetacellToolkit/cli/MCAT options

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

       -g     gamma, graining level of data 
              Proportion of number of single cells in the initial dataset to the number of metacells in the final dataset
              When using MetaCell this correspond to a target gamma (obtained gamma slightly lower)
          
       -s     output, desired metacell file format in output, either 'adata' for a h5ad file or a 'seurat' for a rds file. 
              Output file name will be  'mc_'{output_format}. 
              
       -r     reduction_key (only for SEACells, default "X_pca")

       -y     yaml_file (only for MetaCell2, default None and use default options and gene lists)

### 3.2 Metacell identification on Cd34+ cells using SuperCell

Here we identify metacells from the h5ad file of Cd34+ cells using SuperCell and save the results in a h5ad file. We use 50 principal components, 30 neighbors for the knn and a graining level of 75.

    MCAT -t SuperCell -i  data/cd34_multiome_rna.h5ad -o MCAT_output/SuperCell/cd34/ -n 50 -f 2000 -k 30 -g 75 -s adata

### 3.2 Metacell identification on PBMCs using SEACells

Here we identify metacells from the rds file of PBMCs using SEACells and save the results in a rds file. We use here a graining level of 50.

    MCAT -t SEACells -i data/pbmc.rds -o MCAT_output/SEACells/pbmc/ -n 50 -f 2000 -k 30 -g 50 -s seurat

### 3.3 Metacell identification on PBMCs using MetaCell

Here we identify metacells from the rds file of PBMCs using SEACells and save the results in a h5ad file.

    MCAT -t MetaCell -i data/pbmc.rds -o MCAT_output/MetaCell/pbmc/ -g 50 -s seurat

MetaCell does not use a knn graph from PCA based the highly variable genes but has its own parameters (including different gene list) you can set using a yaml config file and the `-y` argument. You have an example of a such yaml file [here]((/cli/config/MetaCell2_config.yml)) containing the default settings of MetaCell proposed by the authors.

### 3.4 Supervised Metacell identification

You can identify metacells according to a given annotation (e.g. cell types, samples) present in the metadata of the object using the `-a` argument. You can specify a minimum number of metacells to identify (per annotation) using the `-m` argument.

    #Using SEACells
    python cli/SEACellsCL.py -i data/cd34_multiome_rna.h5ad -o testCLI/SEACells_per_celltype_min_5_MC/cd34_multiome_rna/input_raw_adata/ -a celltype -m 5 -n 50 -f 2000 -k 30 -g 75 -s adata

With SuperCell it is possible to use parallel processing using the `-l` argument which gives the number of cores to use.

    #SuperCell parallel metacell identification in each cell type
    Rscript cli/SuperCellCL.R -i data/cd34_multiome_rna.adata -o testCLI/SuperCell_per_celltype/cd34_multiome_rna/input_raw_adata/ -n 50 -f 2000 -k 30 -g 75 -s adata -a celltype -l 6

## Quality control visualization

-   Using scanpy and SEACells (jupyter notebook TO COME).
-   Using Seurat and MetacellToolkit QC (Rmarkdown TO COME)

## Advanced analysis

-   [Analysis](/examples/HLCA_core_atlas.Rmd) of the core HLCA atlas comprising 500'000 cells at the metacell level using SuperCell in command line and STACAS semi-supervised integration.
