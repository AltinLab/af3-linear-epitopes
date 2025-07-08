# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     custom_cell_magics: kql
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: linear-epitope
#     language: python
#     name: python3
# ---

# %%
# %load_ext autoreload
# %autoreload 2


# %%
import polars as pl

# this one is my package
from mdaf3.AF3OutputParser import AF3Output
from mdaf3.FeatureExtraction import *
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import roc_curve, auc
from pathlib import Path
from af3_linear_epitopes import statistics as st
from af3_linear_epitopes import statistics_focal as stf
from scipy import stats


# import dataframes from "staged" directory

fp_in_class_dat = pl.read_parquet(
    "/scratch/sromero/af3-linear-epitopes/data/in_class/focal_protein/staged/in_class_focal_protein.filt.clust.parquet"
).filter(pl.col("representative"))
in_class_dat = pl.read_parquet(
    "/scratch/sromero/af3-linear-epitopes/data/in_class/peptide/staged/in_class_peptide.filt.parquet"
)

all_statistics_in_class_fp = pl.read_parquet(
    "../../data/in_class/focal_protein/staged/01_in_class.exploded.parquet"
)

# %%
all_statistics_in_class_fp = st.pl_structure_fp(
    all_statistics_in_class_fp,
    "../../data/in_class/focal_protein/inference",
)

   # %%
   all_statistics_in_class_fp = all_statistics_in_class_fp.with_columns(
        (pl.col("helix") / 30).alias("helix_percentage"),
        (pl.col("beta") / 30).alias("beta_sheet_percentage"),
        (pl.col("loop") / 30).alias("loop_percentage"),
    )

# %%
all_statistics_in_class_fp

# %%
y_hat_structure_fp = st.normalized_pLDDT_30mer(
    all_statistics_in_class_fp, "helix_percentage", -1
)
y_true = all_statistics_in_class_fp.select(pl.col("epitope"))
in_class_helix_30mer_fp = st.plot_auc_roc_curve(
    y_true,
    y_hat_structure_fp,
    "in_class inverse normalized helix percentage 30-mer fp ROC",
)

# %%
y_hat_structure_fp = st.normalized_pLDDT_30mer(
    all_statistics_in_class_fp, "beta_sheet_percentage", 0
)
y_true = all_statistics_in_class_fp.select(pl.col("epitope"))
in_class_beta_sheet_30mer_fp = st.plot_auc_roc_curve(
    y_true,
    y_hat_structure_fp,
    "in_class inverse normalized beta sheet percentage 30-mer fp ROC",
)

# %%
y_hat_structure_fp = st.normalized_pLDDT_30mer(
    all_statistics_in_class_fp, "loop_percentage", -1
)
y_true = all_statistics_in_class_fp.select(pl.col("epitope"))
in_class_loop_30mer_fp = st.plot_auc_roc_curve(
    y_true, y_hat_structure_fp, "iv_class normalized loop percentage 30-mer fp ROC"
)

# %%
aggregate_mean_pLDDT_fp = all_statistics_in_class_fp.group_by("fp_job_names").agg(
    (pl.col("helix_percentage")).alias("score"),
    pl.col("epitope").alias("epitope"),
)
import sklearn

aggregate_mean_pLDDT_fp = aggregate_mean_pLDDT_fp.with_columns(
    pl.struct(pl.col("score").alias("y_hat"), pl.col("epitope").alias("y_true"))
    .map_elements(lambda x: sklearn.metrics.roc_auc_score(x["y_true"], x["y_hat"]))
    .alias("AUC")
)
mean_auc_pLDDT = aggregate_mean_pLDDT_fp.select("AUC").mean()
print(mean_auc_pLDDT)

# %%
aggregate_mean_pLDDT_fp = all_statistics_in_class_fp.group_by("fp_job_names").agg(
    (pl.col("beta_sheet_percentage")).alias("score"),
    pl.col("epitope").alias("epitope"),
)
import sklearn

aggregate_mean_pLDDT_fp = aggregate_mean_pLDDT_fp.with_columns(
    pl.struct(pl.col("score").alias("y_hat"), pl.col("epitope").alias("y_true"))
    .map_elements(lambda x: sklearn.metrics.roc_auc_score(x["y_true"], x["y_hat"]))
    .alias("AUC")
)
mean_auc_pLDDT = aggregate_mean_pLDDT_fp.select("AUC").mean()
print(mean_auc_pLDDT)

# %%
aggregate_mean_pLDDT_fp = all_statistics_in_class_fp.group_by("fp_job_names").agg(
    (pl.col("loop_percentage")).alias("score"),
    pl.col("epitope").alias("epitope"),
)

aggregate_mean_pLDDT_fp = aggregate_mean_pLDDT_fp.with_columns(
    pl.struct(pl.col("score").alias("y_hat"), pl.col("epitope").alias("y_true"))
    .map_elements(lambda x: sklearn.metrics.roc_auc_score(x["y_true"], x["y_hat"]))
    .alias("AUC")
)
mean_auc_pLDDT = aggregate_mean_pLDDT_fp.select("AUC").mean()
print(mean_auc_pLDDT)

# %%
all_statistics_in_class_fp = st.pl_amino_acids(
    all_statistics_in_class_fp,
    "../../data/in_class/focal_protein/inference",
)

# %%
all_statistics_in_class_fp = all_statistics_in_class_fp.unnest("amino_acid_count")

    # %%
    amino_acid_count_epitope = {
        "A": [],  # Alanine
        "R": [],  # Arginine
        "N": [],  # Asparagine
        "D": [],  # Aspartic Acid
        "C": [],  # Cysteine
        "Q": [],  # Glutamine
        "E": [],  # Glutamic Acid
        "G": [],  # Glycine
        "H": [],  # Histidine
        "I": [],  # Isoleucine
        "L": [],  # Leucine
        "K": [],  # Lysine
        "M": [],  # Methionine
        "F": [],  # Phenylalanine
        "P": [],  # Proline
        "S": [],  # Serine
        "T": [],  # Threonine
        "W": [],  # Tryptophan
        "Y": [],  # Tyrosine
        "V": [],  # Valine
    }
    amino_acid_count_non_epitope = {
        "A": [],  # Alanine
        "R": [],  # Arginine
        "N": [],  # Asparagine
        "D": [],  # Aspartic Acid
        "C": [],  # Cysteine
        "Q": [],  # Glutamine
        "E": [],  # Glutamic Acid
        "G": [],  # Glycine
        "H": [],  # Histidine
        "I": [],  # Isoleucine
        "L": [],  # Leucine
        "K": [],  # Lysine
        "M": [],  # Methionine
        "F": [],  # Phenylalanine
        "P": [],  # Proline
        "S": [],  # Serine
        "T": [],  # Threonine
        "W": [],  # Tryptophan
        "Y": [],  # Tyrosine
        "V": [],  # Valine
    }
    amino_acid_p_values = {
        "A": 0.0,  # Alanine
        "R": 0.0,  # Arginine
        "N": 0.0,  # Asparagine
        "D": 0.0,  # Aspartic Acid
        "C": 0.0,  # Cysteine
        "Q": 0.0,  # Glutamine
        "E": 0.0,  # Glutamic Acid
        "G": 0.0,  # Glycine
        "H": 0.0,  # Histidine
        "I": 0.0,  # Isoleucine
        "L": 0.0,  # Leucine
        "K": 0.0,  # Lysine
        "M": 0.0,  # Methionine
        "F": 0.0,  # Phenylalanine
        "P": 0.0,  # Proline
        "S": 0.0,  # Serine
        "T": 0.0,  # Threonine
        "W": 0.0,  # Tryptophan
        "Y": 0.0,  # Tyrosine
        "V": 0.0,  # Valine
    }
    amino_acid_t_values = {
        "A": 0.0,  # Alanine
        "R": 0.0,  # Arginine
        "N": 0.0,  # Asparagine
        "D": 0.0,  # Aspartic Acid
        "C": 0.0,  # Cysteine
        "Q": 0.0,  # Glutamine
        "E": 0.0,  # Glutamic Acid
        "G": 0.0,  # Glycine
        "H": 0.0,  # Histidine
        "I": 0.0,  # Isoleucine
        "L": 0.0,  # Leucine
        "K": 0.0,  # Lysine
        "M": 0.0,  # Methionine
        "F": 0.0,  # Phenylalanine
        "P": 0.0,  # Proline
        "S": 0.0,  # Serine
        "T": 0.0,  # Threonine
        "W": 0.0,  # Tryptophan
        "Y": 0.0,  # Tyrosine
        "V": 0.0,  # Valine
    }

# %%
filtered_non_epitope = all_statistics_in_class_fp.filter(~pl.col("epitope"))
filtered = all_statistics_in_class_fp.filter(pl.col("epitope"))

symbols = list(amino_acid_count_epitope.keys())
for i in range(0, len(symbols)):
    amino_acid_count_epitope[symbols[i]] = filtered[symbols[i]].to_list()
    amino_acid_count_non_epitope[symbols[i]] = filtered_non_epitope[
        symbols[i]
    ].to_list()
    t_stat, p_value = stats.ttest_ind(
        amino_acid_count_epitope[symbols[i]], amino_acid_count_non_epitope[symbols[i]]
    )
    amino_acid_p_values[symbols[i]] = p_value
    amino_acid_t_values[symbols[i]] = t_stat

# %%
amino_acid_t_values

# %%
st.plot_dictionary_bar_chart(
    amino_acid_p_values,
    "p_values of amino acids tested from epitope to non-epitope regions",
    "Amino Acid Symbols",
    "p-values",
)

# %%
st.plot_dictionary_bar_chart(
    amino_acid_t_values,
    "t_values of amino acids tested from epitope to non-epitope regions",
    "Amino Acid Symbols",
    "t-values",
)
