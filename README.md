# MetacellToolkit

Toolkit for metacell analysis. It consists of a command line tool to easily identify metacells with either SEACells, SuperCell or MetaCell2 (TO COME) with an accompanying R package for metacell quality control and visualization (TO COME).

# Installation

    conda env create -f MetacellToolkit_env.yml
    conda activate MetacellToolkit
    Rscript install.R

# Usage

Command line tool takes as input/output either Seurat .rds or Anndata .h5ad objects.

first we download CD34+ scRNA-seq dataset from Dana's Peer lab

    mkdir data/
    wget https://dp-lab-data-public.s3.amazonaws.com/SEACells-multiome/cd34_multiome_rna.h5ad -O data/cd34_multiome_rna.h5ad

## Metacell identification with SuperCell

    # Print command line arguments
    Rscript cli/SuperCellCL.R -h

    # input raw adata output adata
    Rscript cli/SuperCellCL.R -i data/cd34_multiome_rna.h5ad -o testCLI/SuperCell/cd34_multiome_rna/input_raw_adata/ -n 50 -f 2000 -k 30 -g 75 -s adata

    # input raw adata output seurat
    Rscript cli/SuperCellCL.R -i data/cd34_multiome_rna.h5ad -o testCLI/SuperCell/cd34_multiome_rna/input_raw_adata/ -n 50 -f 2000 -k 30 -g 75 -s seurat

## Metacell identification with SEACells

    # Print command line arguments
    python cli/SEACellsCL.py -h

    # input raw adata output adata
    python cli/SEACellsCL.py -i data/cd34_multiome_rna.h5ad -o test_cli/SEACells/cd34_multiome_rna/ -n 50 -f 2000 -k 30 -g 100 -s adata

    # input raw adata output seurat
    python SEACellsCL.py -i data/cd34_multiome_rna.h5ad -o test_cli/SEACells/cd34_multiome_rna/ -n 50 -f 2000 -k 30 -g 100 -s seurat

## Supervised Metacell identification

You can identifify metacells according to a given annotation (e.g. cell types, samples) present in the metadata of the object using the -a argument. You can specify a minimum number of metacells to identify (per annotation) using the -m argument.

    #Using SEACells
    python cli/SEACellsCL.py -i data/cd34_multiome_rna.h5ad -o testCLI/SEACells_per_celltype_min_5_MC/cd34_multiome_rna/input_raw_adata/ -a celltype -m 5 -n 50 -f 2000 -k 30 -g 75 -s adata

With SuperCell it is possible to use parallel processing for this using the -l argument which gives the number of cores to use

    #SuperCell parallel metacell identification in each cell type
    Rscript cli/SuperCellCL.R -i data/cd34_multiome_rna.adata -o testCLI/SuperCell_per_celltype/cd34_multiome_rna/input_raw_adata/ -n 50 -f 2000 -k 30 -g 75 -s adata -a celltype -l 6

## Main bash CLI
    
    #Allow to execute
    chmod a+x cli/MetacellToolkit.sh
    #Print help 
    cli/MetacellToolkit.sh -h
    # Using SuperCEll
    cli/MetacellToolkit.sh -t SuperCell i data/cd34_multiome_rna.rds - test_output/SuperCell/cd34_multiome_rna/input_raw_rds/ -n 50 -f 2000 -k 30 -g 75 -s adata


## COMING SOON:

-   MetaCell2

# Quality control visualization

-   Using scanpy and SEACells (jupyter notebook TO COME).
-   Using Seurat and MetacellToolkit QC (Rmarkdown TO COME)

# Advanced examples

-   Regulons analysis in CD34+ cells at the metacell level using SCENIC (TO COME)
-   [Analysis](/examples/HLCA_core_atlas.Rmd) of the core HLCA atlas comprising 500'000 cells at the metacell level using SuperCell in command line and STACAS semi-supervised integration.
