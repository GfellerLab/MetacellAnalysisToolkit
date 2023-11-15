#!/bin/bash

echo "creating MetacellAnalysisToolkit conda env"
conda env create -f MetacellAnalysisToolkit_env.yml
conda activate MetacellAnalysisToolkit

echo "Installing additionnal R packages"
Rscript install.R

echo "Making command line script exectuable"
chmod a+x cli/MCAT
chmod a+x cli/SuperCellCL.R 
chmod a+x cli/SEACells.py
chmod a+x cli/MetaCell2CL.py

echo "Installation finished" 
echo "you migh want to add the command line tool path to your PATH environment variable adding the following line to you bashrc" 
