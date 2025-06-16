# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.2
# ---

# %% [markdown]
# ### Overall data structure:
#
# "Focal proteins" are the original proteins where the 30-mer peptides came from
# "Peptides" refer to the 30-mer peptides that either are or are not epitopes (reactive)
#
# Inference (predicted AF3 structures and confidences) is available for both peptides and focal proteins

# %%
# %load_ext autoreload
# %autoreload 2

# %%
import polars as pl

# this one is my package
from mdaf3.AF3OutputParser import AF3Output
from mdaf3.FeatureExtraction import *


from pathlib import Path

# import dataframes from "staged" directory
fp_test_dat = pl.read_parquet(
    "/scratch/sromero/af3-linear-epitopes/data/hv/focal_protein/staged/00_focal_protein.filt.parquet"
)
peptide_test_dat = pl.read_parquet(
    "/scratch/sromero/af3-linear-epitopes/data/hv/peptide/staged/00_hv.filt.parquet"
)

# just hardcoding the path to test data directory
# can swap this out with real data later
DATA_DIR = Path("/scratch/sromero/af3-linear-epitopes/data/test")

# %% [markdown]
# ## Dataframe structure

# %% [markdown]
# ### Focal protein dataframe
#
# Cols: 
# - job_name: name of the af3 job, use it to retreive af3 output object
# - seq: the full protein sequence from which the peptides (30mers) came
# - raw_protein_id: john's original identifier for the protein (stored in list because original data contains some dupes)

# %%
fp_test_dat

# %% [markdown]
# ### Peptide dataframe
#
# Cols: 
# - job_name: name of the af3 job, use it to retreive af3 output object
# - peptide: the 30mer sequence
# - raw_protein_id: john's original identifier for the peptide
# - epitope: whether or not john's assay (pepseq) detected reactivity against the peptide
# - fp_job_names: list of job_names in the focal proteins dataframe that contain this row's 30-mer peptide
# - fp_seq_idxs: at each index in the list fp_job_names, this list contains a list of integer offsets where the 30-mer peptide appears 

# %%
peptide_test_dat

# %%
from af3_linear_epitopes import statistics as st

all_statistics = st.statistics(
    peptide_test_dat,
    "/tgen_labs/altin/alphafold3/runs/linear_peptide/data/hv/peptide/inference",
)


"""
st.bar_graph(
    all_statistics,
    "/scratch/sromero/af3-linear-epitopes/data/test/focal_protein/inference",
)"""

# %%
print(all_statistics)

# %%
true_mean = (
    all_statistics.filter(pl.col("epitope"))
    .select(pl.col("Mean_pLDDT"))
    .to_series()
    .to_list()
)

avg_true_mean = sum(true_mean) / len(true_mean)

false_mean = (
    all_statistics.filter(~pl.col("epitope"))
    .select(pl.col("Mean_pLDDT"))
    .to_series()
    .to_list()
)
avg_false_mean = sum(false_mean) / len(false_mean)

# %%
true_std = (
    all_statistics.filter(pl.col("epitope"))
    .select(pl.col("Std_pLDDT"))
    .to_series()
    .to_list()
)
avg_true_std = sum(true_std) / len(true_std)

false_std = (
    all_statistics.filter(~pl.col("epitope"))
    .select(pl.col("Std_pLDDT"))
    .to_series()
    .to_list()
)
avg_false_std = sum(false_std) / len(false_std)

# %%
true_min = (
    all_statistics.filter(pl.col("epitope"))
    .select(pl.col("Min_pLDDT"))
    .to_series()
    .to_list()
)

avg_true_min = sum(true_min) / len(true_min)

false_min = (
    all_statistics.filter(~pl.col("epitope"))
    .select(pl.col("Min_pLDDT"))
    .to_series()
    .to_list()
)
avg_false_min = sum(false_min) / len(false_min)


print(
    "Average mean,min,std pLLDT values respectively for epitope: "
    + str(avg_true_mean)
    + ", "
    + str(avg_true_min)
    + ", "
    + str(avg_true_std)
    + "\nAverage mean,min,std pLLDT values respectively for non-epitope: "
    + str(avg_false_mean)
    + ", "
    + str(avg_false_min)
    + ", "
    + str(avg_false_std)
)
