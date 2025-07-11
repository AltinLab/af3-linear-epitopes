#!/bin/bash
#SBATCH --job-name=clean_raw_data
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 1
#SBATCH --time=1:00:00
#SBATCH --output=tmp/nextflow/in_class/clean_raw_data.%j.log

. ./scripts/setup.sh

# env vars
export NXF_LOG_FILE=tmp/nextflow/in_class/clean_raw_data/nextflow.log
export NXF_CACHE_DIR=tmp/nextflow/in_class/clean_raw_data/

conda run -n nf-core --live-stream nextflow run \
    ./workflows/00_clean_raw_data.in_class.nf \
        --dset_name in_class