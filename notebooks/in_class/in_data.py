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
import polars as pl

# %%
# this one is my package
from mdaf3.AF3OutputParser import AF3Output
from mdaf3.FeatureExtraction import *
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import roc_curve, auc
from pathlib import Path
from af3_linear_epitopes import statistics as st
from af3_linear_epitopes import statistics_focal as stf


# %% [markdown]
# import dataframes from "staged" directory

# %%
fp_in_class_dat = pl.read_parquet(
    "/scratch/sromero/af3-linear-epitopes/data/in_class/focal_protein/staged/in_class_focal_protein.filt.clust.parquet"
).filter(pl.col("representative"))
in_class_dat = pl.read_parquet(
    "/scratch/sromero/af3-linear-epitopes/data/in_class/peptide/staged/in_class_peptide.filt.parquet"
)

# %%
all_statistics_in_class_fp = pl.read_parquet(
    "../../data/in_class/focal_protein/staged/01_in_class.exploded.parquet"
)


# %%
in_class_dat


# %%
fp_in_class_dat


# %%
all_statistics_in_class_fp


# %%
y = all_statistics_in_class_fp.select("epitope").to_series()
print(len(y))
y_hat_fp_30mer = st.normalized_pLDDT_30mer(
    all_statistics_in_class_fp, "mean_pLDDT_slice"
)
print(len(y_hat_fp_30mer))
in_class_norm_30mer_pLDDT_ROC = st.plot_auc_roc_curve(
    y, y_hat_fp_30mer, "in_class Normalized focal protein pLDDT mean 30-mer ROC"
)
in_class_norm_30mer_pLDDT_ROC.savefig(
    "../../results/figures/in_class_norm_30mer_pLDDT_ROC.png"
)


# %%
fp_aggrigate_30mer = all_statistics_in_class_fp.group_by("peptide").agg(
    (pl.col("mean_pLDDT_slice").mean()).alias("score"),
    pl.col("epitope").first().alias("epitope"),
)


# %%
y_hat_geometric_fp = st.normalized_pLDDT_30mer(fp_aggrigate_30mer, "score")
y_true = fp_aggrigate_30mer.select(pl.col("epitope"))
in_class_geo_mean_9mer_pLDDT_ROC = st.plot_auc_roc_curve(
    y_true, y_hat_geometric_fp, "in_class geometric mean 9-mer fp ROC"
)
in_class_geo_mean_9mer_pLDDT_ROC.savefig(
    "../../results/figures/in_class_geo_mean_9mer_pLDDT_ROC.png"
)


# %%
all_statistics_in_class_fp = all_statistics_in_class_fp.with_columns(
    (pl.col("pLDDT_slice_9mer").list.eval((pl.element() / 100).log().mean().exp()))
    .list.first()
    .alias("Geometric_mean_9mer")
)


# %%
y_hat_geometric = st.normalized_pLDDT_30mer(
    all_statistics_in_class_fp, "Geometric_mean_9mer"
)
in_class_norm_geo_mean_9mer_pLDDT_ROC = st.plot_auc_roc_curve(
    y, y_hat_geometric, "in_class Normalized geometric mean 9-mer fp ROC"
)
in_class_norm_geo_mean_9mer_pLDDT_ROC.savefig(
    "../../results/figures/in_class_norm_geo_mean_9mer_pLDDT_ROC.png"
)


# %%
fp_aggrigate = all_statistics_in_class_fp.group_by("peptide").agg(
    (pl.col("Geometric_mean_9mer").mean()).alias("score"),
    pl.col("epitope").first().alias("epitope"),
)


# %%
y_hat_geometric_fp = fp_aggrigate.select(pl.col("score"))
y_true = fp_aggrigate.select(pl.col("epitope"))
in_class_geo_mean_30mer_pLDDT_ROC = st.plot_auc_roc_curve(
    y_true, y_hat_geometric_fp, "in_class geometric mean 30-mer fp ROC"
)
in_class_geo_mean_30mer_pLDDT_ROC.savefig(
    "../../results/figures/in_class_geo_mean_30mer_pLDDT_ROC.png"
)


# %%
fp_aggrigate_9mer = all_statistics_in_class_fp.group_by("job_name").agg(
    (pl.col("Geometric_mean_9mer")).alias("score"), pl.col("epitope").alias("epitope")
)


# %%
fp_aggrigate_9mer = fp_aggrigate_9mer.with_columns(
    pl.col("score").list.len().alias("#_of_fp_its_in")
)

# %%
percentage = (
    fp_aggrigate_9mer.filter(pl.col("#_of_fp_its_in") > 1).height
    / fp_aggrigate_9mer.height
) * 100
print(percentage)


# %%
import sklearn

# %%
fp_aggrigate_9mer = fp_aggrigate_9mer.with_columns(
    pl.struct(pl.col("score").alias("y_hat"), pl.col("epitope").alias("y_true"))
    .map_elements(lambda x: sklearn.metrics.roc_auc_score(x["y_true"], x["y_hat"]))
    .alias("AUC")
)


# %%
mean_auc = fp_aggrigate_9mer.select("AUC").mean()
print(mean_auc)
