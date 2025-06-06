#!/bin/bash
#SBATCH --job-name=filt_annot
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH -c 1
#SBATCH --time=1:00:00
#SBATCH --output=tmp/nextflow/filt_annot.%j.log

. ./scripts/setup.sh

nextflow run \
    ./workflows/01_filt_annot.nf \
        -resume