#!/bin/bash
#SBATCH --job-name=clean_raw_data
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 1
#SBATCH --time=1:00:00
#SBATCH --output=tmp/nextflow/hv_class/clean_raw_data.%j.log

# env vars
export NXF_LOG_FILE=tmp/nextflow/hv_class/clean_raw_data/nextflow.log
export NXF_CACHE_DIR=tmp/nextflow/hv_class/clean_raw_data/

conda run -n nf-core --live-stream nextflow run \
    ./workflows/00_clean_raw_data.hv_class.nf \
        --dset_name hv_class \
        -resume