#!/bin/bash
#SBATCH --job-name=filt_annot
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 1
#SBATCH --time=1:00:00
#SBATCH --output=tmp/nextflow/hv/filt_annot.%j.log

. ./scripts/setup.sh

# env vars
export NXF_LOG_FILE=tmp/nextflow/hv/filt_annot/nextflow.log
export NXF_CACHE_DIR=tmp/nextflow/hv/filt_annot/


nextflow run \
    ./workflows/01_filt_annot.nf \
        -resume