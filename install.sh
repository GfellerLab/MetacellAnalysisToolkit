conda env create -f MetacellAnalysisToolkit_env.yml
conda activate MetacellAnalysisToolkit
Rscript install.R
chmod a+x cli/MCAT
chmod a+x cli/SuperCellCL.R 
chmod a+x cli/SEACells.py
chmod a+x cli/MetaCell2CL.py