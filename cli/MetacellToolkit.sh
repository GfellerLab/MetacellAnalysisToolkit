#!/bin/bash
set -e

usage() {
  cat <<EOF
usage: $0 options

Constructing metacell from single cell data with SEACells (0.3.3) 'MetaCell2 (0.9.0) or SuperCell (1.0)
Expect a filtered (low quality cells removed) Seurat or Anndata object  

1 - Identifying metacells, 
2 - aggregating counts data per metacell (summing raw counts)
3 - assigning metadata to metacells and computing purities (Assigning metacells to the most aboundant label)

)

OPTIONS:
   -h     Show this message
   
   -t     tool, either 'SEACells', 'MetaCell' or 'SuperCell' 
   
   -i     input_file, either an Anndata object file '.h5ad' or a Seurat object file '.rds' file
   
   -o     outdir, output directory (default ./)

   -n     dims, number of principal components to use (only for SEACells and SuperCell, default 50) 
   
   -f     n_features, number of highly variable genes use to compute the initial PCA (only for SEACells and SuperCell, default 2000) 
   
   -k     k_knn, number of neighbors to construct the knn graph (only for SEACells and SuperCell, default 30)
   
   -g     gamma, 	graining level of data 
          Proportion of number of single cells in the initial dataset to the number of metacells in the final dataset
          When using MetaCell this correspond to a target gamma (obtained gamma slightly lower)
          
   -s     output, desired file format in output for metacells, either an Anndata object file '.h5ad' or a Seurat object file '.rds' file
   
   -r     reduction_key (only for SEACells, default "X_pca")
   
   -y     yaml_file (only for MetaCell2, default None and use default options and gene lists)
EOF
}

SCRIPT_DIR=$(dirname $0)
TOOL_CMD='SuperCell'
OUTDIR==$(dirname $0)
QUIET=""
PREPRO=""
REDUCTION_KEY=""
YAML_FILE=""

while getopts "hr:t:i:o:pn:f:k:g:s:r:y:" OPTION; do
  case $OPTION in
  h)
    usage
    exit 1
    ;;
  t)
    TOOL=$OPTARG
    ;;
  i)
    INPUT_FILE=$OPTARG
    ;;
  o)
    OUTDIR=$OPTARG
    ;;
  n)
    DIMS=$OPTARG
    ;;
  f)
    N_FEATURES=$OPTARG
    ;;
  k)
    K_KNN=$OPTARG
    ;;
  g)
    GAMMA=$OPTARG
    ;;
  s)
    OUTPUT=$OPTARG
    ;;
  r)
    REDUCTION_KEY="-r ${OPTARG}"
    ;;
  y)
    YAML_FILE="-y ${OPTARG}"
    ;;
  ?)
    usage
    exit
    ;;
  esac
done

echo ${TOOL}

echo ${INPUT_FILE}

# Check tool version
if [[ -z $TOOL ]] || (([ "$TOOL" != "SuperCell" ] && [ "$TOOL" != "SEACells" ]) &&  [ "$TOOL" != "MetaCell" ]); then
  echo "Invalid tool $TOOL"
  usage
  exit 1
fi

echo "Identifying metacells..."
if [ "$TOOL" == "SEACells" ]; then
  python  "$SCRIPT_DIR/SEACellsCL.py" -i ${INPUT_FILE} -o ${OUTDIR} -n ${DIMS} -f ${N_FEATURES} -g ${GAMMA} -s ${OUTPUT} ${REDUCTION_KEY}
elif [ "$TOOL" == "MetaCell" ]; then
  python  "$SCRIPT_DIR/MetaCell2CL.py" -i ${INPUT_FILE} -o ${OUTDIR} -g ${GAMMA} -s ${OUTPUT} ${YAML_FILE}
elif [ "$TOOL" == "SuperCell" ]; then
  Rscript "$SCRIPT_DIR/SuperCellCL.R" -i ${INPUT_FILE} -o ${OUTDIR} -n ${DIMS} -f ${N_FEATURES} -g ${GAMMA} -s ${OUTPUT} 
fi

echo "Done."