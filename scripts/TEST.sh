#!/bin/bash
#SBATCH --job-name=test
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 1
#SBATCH --time=1:00:00
#SBATCH --output=tmp/nextflow/test.%j.log

. ./scripts/setup.sh

# # env vars
# export NXF_LOG_FILE=tmp/nextflow/test/clean_raw_data/nextflow.log
# export NXF_CACHE_DIR=tmp/nextflow/test/clean_raw_data/

# nextflow run \
#     ./workflows/00_clean_raw_data.nf \
#         --data_dir ./data/test \
#         -resume

# # env vars
# export NXF_LOG_FILE=tmp/nextflow/test/filt_annot/nextflow.log
# export NXF_CACHE_DIR=tmp/nextflow/test/filt_annot/

# nextflow run \
#     ./workflows/01_filt_annot.nf \
#         --data_dir ./data/test \
#         -resume

# env vars
export NXF_LOG_FILE=tmp/nextflow/test/inference/nextflow.log
export NXF_CACHE_DIR=tmp/nextflow/test/inference/


# nextflow run \
#     ./workflows/02_msa_peptide.nf \
#     --data_dir ./data/test \
#     -resume && \
# nextflow run \
#     ./workflows/02_msa_focal_protein.nf \
#     --data_dir ./data/test \
#     -resume


nextflow run \
    ./workflows/02_inference_peptide.nf \
        --data_dir ./data/test \
        --seeds 1,2 

# nextflow run \
#     ./workflows/02_inference_focal_protein.nf \
#         --data_dir ./data/test \
#          --seeds 1,2 \
#         -resume 