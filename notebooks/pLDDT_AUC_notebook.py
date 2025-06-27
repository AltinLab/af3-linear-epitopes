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
fp_test_dat = pl.read_parquet(
    "/scratch/sromero/af3-linear-epitopes/data/hv/focal_protein/staged/00_focal_protein.filt.parquet"
)
all_statistics = pl.read_parquet(
    "/scratch/sromero/af3-linear-epitopes/data/hv/peptide/staged/01_hv.features.parquet"
)
all_statistics_fp = pl.read_parquet(
    "/scratch/sromero/af3-linear-epitopes/data/hv/peptide/staged/01_hv.exploded.parquet"
)


y_hat = st.normalized_pLDDT_30mer(all_statistics, "Mean_pLDDT")


y = all_statistics.select("epitope").to_series()


st.plot_auc_roc_curve(y, y_hat, "Normalized Mean pLDDT 30-mer ROC")


all_statistics = all_statistics.with_columns(
    pl.col("9mer_Mean_pLDDT").list.min().alias("Min_of_means_9mer_peptide")
)


all_statistics = all_statistics.with_columns(
    pl.col("9mer_Mean_pLDDT").list.max().alias("Max_of_means_9mer_peptide")
)


y_hat_min = st.normalized_pLDDT_30mer(all_statistics, "Min_of_means_9mer_peptide")
st.plot_auc_roc_curve(y, y_hat_min, "Normalized Min pLDDT 9mer ROC")


y_hat_max = st.normalized_pLDDT_30mer(all_statistics, "Max_of_means_9mer_peptide")
st.plot_auc_roc_curve(y, y_hat_max, "Normalized Max pLDDT 9mer ROC")


y_hat_weight = st.normalized_pLDDT_30mer(all_statistics, "atomic_weight")
st.plot_auc_roc_curve(y, y_hat_weight, "Normalized atomic weight 30-mer ROC")


fp_test_dat


all_statistics_fp


y = all_statistics_fp.select("epitope").to_series()
print(len(y))
y_hat_fp_30mer = st.normalized_pLDDT_30mer(all_statistics_fp, "mean_pLDDT_slice")
print(len(y_hat_fp_30mer))
st.plot_auc_roc_curve(
    y, y_hat_fp_30mer, "Normalized focal protein pLDDT mean 30-mer ROC"
)


fp_aggrigate_30mer = all_statistics_fp.group_by("peptide").agg(
    (pl.col("mean_pLDDT_slice").mean()).alias("score"),
    pl.col("epitope").first().alias("epitope"),
)


y_hat_geometric_fp = st.normalized_pLDDT_30mer(fp_aggrigate_30mer, "score")
y_true = fp_aggrigate_30mer.select(pl.col("epitope"))
st.plot_auc_roc_curve(y_true, y_hat_geometric_fp, "geometric mean 9-mer fp ROC")


all_statistics_fp


all_statistics_fp = all_statistics_fp.with_columns(
    (pl.col("pLDDT_slice_9mer").list.eval((pl.element() / 100).log().mean().exp()))
    .list.first()
    .alias("Geometric_mean_9mer")
)


all_statistics_fp


y_hat_geometric = st.normalized_pLDDT_30mer(all_statistics_fp, "Geometric_mean_9mer")
st.plot_auc_roc_curve(y, y_hat_geometric, "Normalized geometric mean 9-mer fp ROC")


fp_aggrigate = all_statistics_fp.group_by("peptide").agg(
    (pl.col("Geometric_mean_9mer").mean()).alias("score"),
    pl.col("epitope").first().alias("epitope"),
)


y_hat_geometric_fp = fp_aggrigate.select(pl.col("score"))
y_true = fp_aggrigate.select(pl.col("epitope"))
st.plot_auc_roc_curve(
    y_true, y_hat_geometric_fp, "Normalized geometric mean 9-mer fp ROC"
)


fp_aggrigate_9mer = all_statistics_fp.group_by("fp_job_names").agg(
    (pl.col("Geometric_mean_9mer")).alias("score"), pl.col("epitope").alias("epitope")
)


fp_aggrigate_9mer


import sklearn

fp_aggrigate_9mer = fp_aggrigate_9mer.with_columns(
    pl.struct(pl.col("score").alias("y_hat"), pl.col("epitope").alias("y_true"))
    .map_elements(lambda x: sklearn.metrics.roc_auc_score(x["y_true"], x["y_hat"]))
    .alias("AUC")
)


mean_auc = fp_aggrigate_9mer.select("AUC").mean()
print(mean_auc)
