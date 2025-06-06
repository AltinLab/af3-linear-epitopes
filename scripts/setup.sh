#!/bin/bash

# load module for nextflow
. /tgen_labs/altin/miniconda3/etc/profile.d/conda.sh

conda activate nf-core

# env vars
export NXF_LOG_FILE=tmp/nextflow/nextflow.log
export NXF_CACHE_DIR=tmp/nextflow