#!/bin/bash
#SBATCH --job-name=extract_conf_iedb_bp3
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --time=5-00:00:00
#SBATCH -c 16
#SBATCH --output=tmp/nextflow/iedb_bp3/focal_protein/extract_conf.%j.log

# env vars
export NXF_LOG_FILE=tmp/nextflow/iedb_bp3/focal_protein/bp3/nextflow.log

conda run -n nf-core --live-stream nextflow run \
    ./workflows/06_extract_conf.nf \
        --dset_name iedb_bp3