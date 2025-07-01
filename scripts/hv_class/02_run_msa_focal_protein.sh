#!/bin/bash
#SBATCH --job-name=msa_focal_protein
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --time=5-00:00:00
#SBATCH -c 16
#SBATCH --output=tmp/nextflow/hv_class/focal_protein/msa.%j.log


# env vars
export NXF_LOG_FILE=tmp/nextflow/hv_class/focal_protein/msa/nextflow.log
export NXF_CACHE_DIR=tmp/nextflow/hv_class/focal_protein/msa/

conda run -n nf-core --live-stream nextflow run \
    ./workflows/02_msa_focal_protein.nf \
        --dset_name hv_class
