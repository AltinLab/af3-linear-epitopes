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


# import dataframes from "staged" directory

fp_hv_class_dat = pl.read_parquet(
    "../../data/hv_class/focal_protein/staged/hv_class_focal_protein.filt.clust.parquet"
).filter(pl.col("representative"))
hv_class_dat = pl.read_parquet(
    "../../data/hv_class/peptide/staged/hv_class_peptide.filt.parquet"
)
all_statistics_hv_class = pl.read_parquet(
    "../../data/hv_class/peptide/staged/01_hv_class.features.parquet"
)
all_statistics_hv_class_fp = pl.read_parquet(
    "../../data/hv_class/focal_protein/staged/01_hv_class.exploded.parquet"
)

# %%
all_statistics_hv_class_fp = st.pl_structure(
    all_statistics_hv_class_fp,
    "../../data/hv_class/peptide/inference",
)

# %%
y_hat_structure_fp = st.normalized_pLDDT_30mer(
    all_statistics_hv_class_fp, "helix_percentage"
)
y_true = all_statistics_hv_class_fp.select(pl.col("epitope"))
hv_class_geometric_mean_pLDDT_30mer_fp = st.plot_auc_roc_curve(
    y_true, y_hat_structure_fp, "hv_class normalized helix percentage 30-mer fp ROC"
)

# %%
y_hat_structure_fp = st.normalized_pLDDT_30mer(
    all_statistics_hv_class_fp, "beta_sheet_percentage"
)
y_true = all_statistics_hv_class_fp.select(pl.col("epitope"))
hv_class_geometric_mean_pLDDT_30mer_fp = st.plot_auc_roc_curve(
    y_true,
    y_hat_structure_fp,
    "hv_class normalized beta sheet percentage 30-mer fp ROC",
)

# %%
y_hat_structure_fp = st.normalized_pLDDT_30mer(
    all_statistics_hv_class_fp, "loop_percentage"
)
y_true = all_statistics_hv_class_fp.select(pl.col("epitope"))
hv_class_geometric_mean_pLDDT_30mer_fp = st.plot_auc_roc_curve(
    y_true, y_hat_structure_fp, "hv_class normalized loop percentage 30-mer fp ROC"
)

# %%
aggregate_mean_pLDDT_fp = all_statistics_hv_class_fp.group_by("fp_job_names").agg(
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
aggregate_mean_pLDDT_fp = all_statistics_hv_class_fp.group_by("fp_job_names").agg(
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
aggregate_mean_pLDDT_fp = all_statistics_hv_class_fp.group_by("fp_job_names").agg(
    (pl.col("loop_percentage")).alias("score"),
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
