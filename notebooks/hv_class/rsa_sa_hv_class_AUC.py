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


# %%
# import dataframes from "staged" directory
fp_hv_class_dat = pl.read_parquet(
    "../../data/hv_class/focal_protein/staged/hv_class_focal_protein.filt.clust.parquet"
).filter(pl.col("representative"))
hv_class_dat = pl.read_parquet(
    "../../data/hv_class/peptide/staged/hv_class_peptide.filt.parquet"
)
all_statistics_hv_class_fp = pl.read_parquet(
    "../../data/hv_class/focal_protein/staged/01_hv_class.exploded.parquet"
)
rsa_data = pl.read_parquet(
    "../../data/hv_class/focal_protein/staged/05_focal_protein.rsa.parquet"
).filter(pl.col("representative"))


# %%
all_statistics_hv_class_fp = rsa_data.join(
    all_statistics_hv_class_fp, left_on="job_name", right_on="fp_job_names"
)


# %%
all_statistics_hv_class_fp = all_statistics_hv_class_fp.drop(
    "seq_right", "raw_protein_ids_right"
)


# %%
all_statistics_hv_class_fp = all_statistics_hv_class_fp.rename(
    {"job_name": "fp_job_name", "job_name_right": "job_name"}
)


# %%
all_statistics_hv_class_fp = st.rsa_mean(all_statistics_hv_class_fp)


# %%
y_hat_RSA_fp = st.normalized_pLDDT_30mer(
    all_statistics_hv_class_fp, "mean_rsa_slice", -1
)
y_true_RSA = all_statistics_hv_class_fp.select(pl.col("epitope"))
inverse_norm_rsa_mean_30mer_fp_hv_class = st.plot_auc_roc_curve(
    y_true_RSA,
    y_hat_RSA_fp,
    "hv_class inverse Normalized mean RSA values for 30mer fp ROC",
)
inverse_norm_rsa_mean_30mer_fp_hv_class.savefig(
    "../../results/figures/inverse_norm_rsa_mean_30mer_fp_hv_class.png"
)


# %%
y_hat_SA_fp = st.normalized_pLDDT_30mer(all_statistics_hv_class_fp, "mean_sa_slice", -1)
y_true_SA = all_statistics_hv_class_fp.select(pl.col("epitope"))
inverse_norm_sa_mean_30mer_fp_hv_class = st.plot_auc_roc_curve(
    y_true_SA,
    y_hat_SA_fp,
    "hv_class inverse Normalized mean SA values for 30mer fp ROC",
)
inverse_norm_sa_mean_30mer_fp_hv_class.savefig(
    "../../results/figures/inverse_norm_sa_mean_30mer_fp_hv_class.png"
)


# %%
aggregate_mean_rsa = all_statistics_hv_class_fp.group_by("peptide").agg(
    (pl.col("mean_rsa_slice").mean()).alias("score"),
    pl.col("epitope").first().alias("epitope"),
)


# %%
y_hat_RSA_fp = aggregate_mean_rsa.select(pl.col("score"))
y_true_RSA = aggregate_mean_rsa.select(pl.col("epitope"))
rsa_mean_aggregate_30mer_fp_hv_class = st.plot_auc_roc_curve(
    y_true_RSA, y_hat_RSA_fp, "hv_class aggregate mean RSA values for 30mer fp ROC"
)
rsa_mean_aggregate_30mer_fp_hv_class.savefig(
    "../../results/figures/rsa_mean_aggregate_30mer_fp_hv_class.png"
)


# %%
aggregate_mean_sa = all_statistics_hv_class_fp.group_by("peptide").agg(
    (pl.col("mean_sa_slice").mean()).alias("score"),
    pl.col("epitope").first().alias("epitope"),
)


# %%
y_hat_SA_fp = aggregate_mean_sa.select(pl.col("score"))
y_true_SA = aggregate_mean_sa.select(pl.col("epitope"))
sa_mean_aggregate_30mer_fp_hv_class = st.plot_auc_roc_curve(
    y_true_SA,
    y_hat_SA_fp,
    "hv_class aggregate mean SA values for 30mer fp ROC",
)
sa_mean_aggregate_30mer_fp_hv_class.savefig(
    "../../results/figures/sa_mean_aggregate_30mer_fp_hv_class.png"
)


# %%
aggregate_mean_rsa_fp = all_statistics_hv_class_fp.group_by("fp_job_name").agg(
    (pl.col("mean_rsa_slice")).alias("score"),
    pl.col("epitope").alias("epitope"),
)
aggregate_mean_sa_fp = all_statistics_hv_class_fp.group_by("fp_job_name").agg(
    (pl.col("mean_sa_slice")).alias("score"),
    pl.col("epitope").alias("epitope"),
)


# %%
import sklearn


# %%
aggregate_mean_rsa_fp = aggregate_mean_rsa_fp.with_columns(
    pl.struct(pl.col("score").alias("y_hat"), pl.col("epitope").alias("y_true"))
    .map_elements(lambda x: sklearn.metrics.roc_auc_score(x["y_true"], x["y_hat"]))
    .alias("AUC")
)
aggregate_mean_sa_fp = aggregate_mean_sa_fp.with_columns(
    pl.struct(pl.col("score").alias("y_hat"), pl.col("epitope").alias("y_true"))
    .map_elements(lambda x: sklearn.metrics.roc_auc_score(x["y_true"], x["y_hat"]))
    .alias("AUC")
)


# %%
mean_auc_rsa = aggregate_mean_rsa_fp.select("AUC").mean()
mean_auc_sa = aggregate_mean_sa_fp.select("AUC").mean()
print(mean_auc_rsa)
print(mean_auc_sa)

# %%
auc_rsa_data = aggregate_mean_rsa_fp["AUC"].to_list()
auc_sa_data = aggregate_mean_sa_fp["AUC"].to_list()
st.display_boxplot(auc_rsa_data, "Mean auc for RSA values box plot", "", "AUC")
st.display_boxplot(auc_sa_data, "Mean auc for SA values box plot", "", "AUC")

# %%
aggregate_mean_sa_fp.filter(pl.col("AUC") == 0)

# %%
zero_auc = aggregate_mean_sa_fp.filter(pl.col("AUC") == 0)
zero_auc.filter(pl.col("score").list.len() == 2)
