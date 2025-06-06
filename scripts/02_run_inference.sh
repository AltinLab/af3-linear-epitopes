#!/bin/bash
#SBATCH --job-name=alphafold
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --time=5-00:00:00
#SBATCH --output=run_alpha.%j.log

module load Java/11.0.2

nextflow run \
    -w /scratch/lwoods/work \
    -c /home/lwoods/workspace/af3-nf/nextflow.config \
    /home/lwoods/workspace/af3-nf/af3_single_prot_msa.nf \
        --input_csv '../jobfiles/HV1_annot_tax_af3job.csv' \
        --seq_col 'peptide' \
        --msa_subdir 'peptide' \
        -resume \
        --msa_db 'https://pub-vscratch.vast.rc.tgen.org' && \
nextflow run \
    -w /scratch/lwoods/work \
    -c /home/lwoods/workspace/af3-nf/nextflow.config \
    /home/lwoods/workspace/af3-nf/af3_single_prot_inference.nf \
        --input_csv '../jobfiles/HV1_annot_tax_af3job.csv' \
        --out_dir '/tgen_labs/altin/alphafold3/runs/linear_peptide/data/1_run_inference' \
        --msa_db 'https://pub-vscratch.vast.rc.tgen.org' \
        --seeds 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100 \
        --compress \
        --check_inf_exists \
        --seq_col 'peptide' \
        --msa_subdir 'peptide' \
        --job_name_col 'job_name' \
        -resume