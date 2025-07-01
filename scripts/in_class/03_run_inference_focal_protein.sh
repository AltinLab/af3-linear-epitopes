#!/bin/bash
#SBATCH --job-name=inference_focal_protein
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --time=5-00:00:00
#SBATCH -c 16
#SBATCH --output=tmp/nextflow/in_class/focal_protein/in_classference.%j.log

. ./scripts/setup.sh

# env vars
export NXF_LOG_FILE=tmp/nextflow/in_class/focal_protein/inference/nextflow.log
export NXF_CACHE_DIR=tmp/nextflow/in_class/focal_protein/inference/

conda run -n nf-core --live-stream nextflow run \
    ./workflows/03_inference_focal_protein.nf \
        --dset_name in_class \
        --seeds 1,2,3,4,5,6,7,8,9,10
