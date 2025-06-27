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
    "../data/hv/focal_protein/staged/00_focal_protein.filt.parquet"
)
all_statistics = pl.read_parquet("../data/hv/peptide/staged/01_hv.features.parquet")
all_statistics_fp = pl.read_parquet("../data/hv/peptide/staged/01_hv.exploded.parquet")


y_hat_weight = st.normalized_pLDDT_30mer(all_statistics, "atomic_weight")
y = all_statistics.select("epitope").to_series()
st.plot_auc_roc_curve(y, y_hat_weight, "Normalized atomic weight 30-mer ROC")


y_hat_mass_fp = st.normalized_pLDDT_30mer(all_statistics, "atomic_weight")
y_true = all_statistics.select(pl.col("epitope"))
st.plot_auc_roc_curve(y_true, y_hat_mass_fp, "Normalized atomic mass for 30-mer fp ROC")


aggrigate_atomic_weights = all_statistics.group_by("peptide").agg(
    (pl.col("9mer_weight").list.max()).first().alias("score"),
    pl.col("epitope").first().alias("epitope"),
)


aggrigate_atomic_weights


y_hat_mass_fp = st.normalized_pLDDT_30mer(aggrigate_atomic_weights, "score")
y_true = aggrigate_atomic_weights.select(pl.col("epitope"))
st.plot_auc_roc_curve(
    y_true, y_hat_mass_fp, "Normalized max of atomic mass for 9-mers ROC"
)


aggrigate_atomic_weights_mean = all_statistics_fp.group_by("fp_job_names").agg(
    (pl.col("9mer_weight").list.mean()).alias("score"),
    pl.col("epitope").alias("epitope"),
)
aggrigate_atomic_weights_min = all_statistics_fp.group_by("fp_job_names").agg(
    (pl.col("9mer_weight").list.min()).alias("score"),
    pl.col("epitope").alias("epitope"),
)
aggrigate_atomic_weights_max = all_statistics_fp.group_by("fp_job_names").agg(
    (pl.col("9mer_weight").list.max()).alias("score"),
    pl.col("epitope").alias("epitope"),
)


import sklearn

aggrigate_atomic_weights_mean = aggrigate_atomic_weights_mean.with_columns(
    pl.struct(pl.col("score").alias("y_hat"), pl.col("epitope").alias("y_true"))
    .map_elements(lambda x: sklearn.metrics.roc_auc_score(x["y_true"], x["y_hat"]))
    .alias("AUC")
)
aggrigate_atomic_weights_min = aggrigate_atomic_weights_min.with_columns(
    pl.struct(pl.col("score").alias("y_hat"), pl.col("epitope").alias("y_true"))
    .map_elements(lambda x: sklearn.metrics.roc_auc_score(x["y_true"], x["y_hat"]))
    .alias("AUC")
)
aggrigate_atomic_weights_max = aggrigate_atomic_weights_max.with_columns(
    pl.struct(pl.col("score").alias("y_hat"), pl.col("epitope").alias("y_true"))
    .map_elements(lambda x: sklearn.metrics.roc_auc_score(x["y_true"], x["y_hat"]))
    .alias("AUC")
)


mean_auc_fp = aggrigate_atomic_weights_mean.select("AUC").mean()
mean_auc_min = aggrigate_atomic_weights_min.select("AUC").mean()
mean_auc_max = aggrigate_atomic_weights_max.select("AUC").mean()
print(mean_auc_fp)
print(mean_auc_min)
print(mean_auc_max)
