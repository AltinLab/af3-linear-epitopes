#!/bin/bash
#SBATCH --job-name=msa_focal_protein
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --time=5-00:00:00
#SBATCH -c 16
#SBATCH --output=tmp/nextflow/in_class/focal_protein/msa.%j.log

. ./scripts/setup.sh

# env vars
export NXF_LOG_FILE=tmp/nextflow/in_class/focal_protein/msa/nextflow.log
export NXF_CACHE_DIR=tmp/nextflow/in_class/focal_protein/msa/

conda run -n nf-core --live-stream nextflow run \
    ./workflows/02_msa_focal_protein.nf \
        --dset_name in_class
