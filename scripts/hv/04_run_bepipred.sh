#!/bin/bash
#SBATCH --job-name=bp3_focal_protein
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --time=5-00:00:00
#SBATCH -c 16
#SBATCH --output=tmp/nextflow/hv/focal_protein/bp3.%j.log

. ./scripts/setup.sh

# env vars
export NXF_LOG_FILE=tmp/nextflow/hv/focal_protein/bp3/nextflow.log

nextflow run \
    ./workflows/04_bepipred_focal_protein.nf \
        --dset_name hv 