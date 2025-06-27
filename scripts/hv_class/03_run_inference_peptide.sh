#!/bin/bash
#SBATCH --job-name=inference_peptide
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --time=5-00:00:00
#SBATCH -c 16
#SBATCH --output=tmp/nextflow/hv_class/peptide/inference.%j.log

. ./scripts/setup.sh

# env vars
export NXF_LOG_FILE=tmp/nextflow/hv_class/peptide/inference/nextflow.log
export NXF_CACHE_DIR=tmp/nextflow/hv_class/peptide/inference/

nextflow run \
    ./workflows/03_inference_peptide.nf \
        --dset_name hv_class \
        --seeds 1,2,3,4,5,6,7,8,9,10 
