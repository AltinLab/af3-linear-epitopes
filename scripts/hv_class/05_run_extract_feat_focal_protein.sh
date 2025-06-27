#!/bin/bash
#SBATCH --job-name=extract_feat
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --time=5-00:00:00
#SBATCH -c 16
#SBATCH --output=tmp/nextflow/hv/focal_protein/extract_feat.%j.log

. ./scripts/setup.sh

# env vars
export NXF_LOG_FILE=tmp/nextflow/hv/focal_protein/extract_feat/nextflow.log
export NXF_CACHE_DIR=tmp/nextflow/hv/focal_protein/extract_feat/

nextflow run \
    ./workflows/05_extract_feat_focal_protein.nf \
        --dset_name hv
