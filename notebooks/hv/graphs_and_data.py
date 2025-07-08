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


pLDDT_statistics_30mer = st.plot_epitope_non_epitope_stats_30mer(
    avg_true_mean,
    avg_true_min,
    avg_true_std,
    avg_false_mean,
    avg_false_min,
    avg_false_std,
)


true_mean_min_9mer = (
    all_statistics.filter(pl.col("epitope"))
    .select(pl.col("9mer_Mean_pLDDT").list.min())
    .to_series()
    .to_list()
)

avg_true_mean_min_9mer = sum(true_mean_min_9mer) / len(true_mean_min_9mer)

false_mean_min_9mer = (
    all_statistics.filter(~pl.col("epitope"))
    .select(pl.col("9mer_Mean_pLDDT").list.min())
    .to_series()
    .to_list()
)
avg_false_mean_min_9mer = sum(false_mean_min_9mer) / len(false_mean_min_9mer)

# --------------------------------------------------------------------------------------
true_min_min_9mer = (
    all_statistics.filter(pl.col("epitope"))
    .select(pl.col("9mer_min_pLDDT").list.min())
    .to_series()
    .to_list()
)

avg_true_min_min_9mer = sum(true_min_min_9mer) / len(true_min_min_9mer)

false_min_min_9mer = (
    all_statistics.filter(~pl.col("epitope"))
    .select(pl.col("9mer_min_pLDDT").list.min())
    .to_series()
    .to_list()
)
avg_false_min_min_9mer = sum(false_min_min_9mer) / len(false_min_min_9mer)

# -------------------------------------------------------------------------------------
true_std_min_9mer = (
    all_statistics.filter(pl.col("epitope"))
    .select(pl.col("9mer_std_pLDDT").list.min())
    .to_series()
    .to_list()
)

avg_true_std_min_9mer = sum(true_std_min_9mer) / len(true_std_min_9mer)

false_std_min_9mer = (
    all_statistics.filter(~pl.col("epitope"))
    .select(pl.col("9mer_std_pLDDT").list.min())
    .to_series()
    .to_list()
)
avg_false_std_min_9mer = sum(false_std_min_9mer) / len(false_std_min_9mer)


print(
    "Average mean,min,std pLLDT values respectively for epitope: "
    + str(avg_true_mean_min_9mer)
    + ", "
    + str(avg_true_min_min_9mer)
    + ", "
    + str(avg_true_std_min_9mer)
    + "\nAverage mean,min,std pLLDT values respectively for non-epitope: "
    + str(avg_false_mean_min_9mer)
    + ", "
    + str(avg_false_min_min_9mer)
    + ", "
    + str(avg_false_std_min_9mer)
)


pLDDT_avg_9mer = st.plot_epitope_non_epitope_stats_9mer(
    avg_true_mean_min_9mer,
    avg_true_std_min_9mer,
    avg_false_mean_min_9mer,
    avg_false_std_min_9mer,
)


y_hat = st.normalized_pLDDT_30mer(all_statistics, "Mean_pLDDT")


y = all_statistics.select("epitope").to_series()


mean_of_30mer = st.plot_auc_roc_curve(y, y_hat, "Normalized Mean of pLDDT 30-mear ROC")


all_statistics = all_statistics.with_columns(
    pl.col("9mer_Mean_pLDDT").list.min().alias("Min_of_means_9mer_peptide")
)


all_statistics = all_statistics.with_columns(
    pl.col("9mer_Mean_pLDDT").list.max().alias("Max_of_means_9mer_peptide")
)


y_hat_min = st.normalized_pLDDT_30mer(all_statistics, "Min_of_means_9mer_peptide")
min_roc_9mer = st.plot_auc_roc_curve(y, y_hat_min, "Normalized Min of pLDDT 9-mer ROC")


y_hat_max = st.normalized_pLDDT_30mer(all_statistics, "Max_of_means_9mer_peptide")
max_roc_9mer = st.plot_auc_roc_curve(y, y_hat_max, "Normalized Max of pLDDT 9-mer ROC")


y_hat_weight = st.normalized_pLDDT_30mer(all_statistics, "atomic_weight")
norm_atomic_weight_30mer = st.plot_auc_roc_curve(
    y, y_hat_weight, "Normalized atomic weight 30-mer peptide ROC"
)


y = all_statistics_fp.select("epitope").to_series()
y_hat_fp_30mer = st.normalized_pLDDT_30mer(all_statistics_fp, "mean_pLDDT_slice")
norm_fp_pLDDT_mean = st.plot_auc_roc_curve(
    y, y_hat_fp_30mer, "Normalized focal protein pLDDT mean 30-mer ROC"
)


# takes the mean of every mean pLDDT score for a peptide in different focal protiens
fp_aggrigate_30mer = all_statistics_fp.group_by("peptide").agg(
    (pl.col("pLDDT_slice_9mer").list.mean().mean()).alias("score"),
    pl.col("epitope").first().alias("epitope"),
)


y_hat_geometric_fp = st.normalized_pLDDT_30mer(fp_aggrigate_30mer, "score")
y_true = fp_aggrigate_30mer.select(pl.col("epitope"))
norm_mean_pLDDT_9mer_fp = st.plot_auc_roc_curve(
    y_true,
    y_hat_geometric_fp,
    "nomralized mean of every mean pLDDT score for a peptide in different focal protiens ROC",
)


all_statistics_fp = all_statistics_fp.with_columns(
    (pl.col("pLDDT_slice_9mer").list.eval((pl.element() / 100).log().mean().exp()))
    .list.first()
    .alias("Geometric_mean_9mer")
)


y_hat_geometric = st.normalized_pLDDT_30mer(all_statistics_fp, "Geometric_mean_9mer")
norm_geometric_mean_9mer_fp = st.plot_auc_roc_curve(
    y, y_hat_geometric, "Normalized geometric mean 9-mer focal protein ROC"
)


fp_aggrigate = all_statistics_fp.group_by("peptide").agg(
    (pl.col("Geometric_mean_9mer").mean()).alias("score"),
    pl.col("epitope").first().alias("epitope"),
)


y_hat_geometric_fp = fp_aggrigate.select(pl.col("score"))
y_true = fp_aggrigate.select(pl.col("epitope"))
geometric_mean_9mer_fp = st.plot_auc_roc_curve(
    y_true, y_hat_geometric_fp, "geometric mean 9-mer fp ROC"
)


fp_aggrigate_9mer = all_statistics_fp.group_by("fp_job_names").agg(
    (pl.col("Geometric_mean_9mer")).alias("score"), pl.col("epitope").alias("epitope")
)


import sklearn

fp_aggrigate_9mer = fp_aggrigate_9mer.with_columns(
    pl.struct(pl.col("score").alias("y_hat"), pl.col("epitope").alias("y_true"))
    .map_elements(lambda x: sklearn.metrics.roc_auc_score(x["y_true"], x["y_hat"]))
    .alias("AUC")
)


mean_auc = fp_aggrigate_9mer.select("AUC").mean()
print(mean_auc)


y_hat_mass_fp = st.normalized_pLDDT_30mer(all_statistics_fp, "atomic_weight")
y_true = all_statistics_fp.select(pl.col("epitope"))
norm_weight_30mer_fp = st.plot_auc_roc_curve(
    y_true, y_hat_mass_fp, "Normalized atomic weight for 30-mer fp ROC"
)


aggrigate_atomic_weights = all_statistics.group_by("peptide").agg(
    (pl.col("9mer_weight").list.min()).first().alias("score"),
    pl.col("epitope").first().alias("epitope"),
)
y_hat_mass_fp = st.normalized_pLDDT_30mer(aggrigate_atomic_weights, "score")
y_true = aggrigate_atomic_weights.select(pl.col("epitope"))
norm_min_weight_9mer_fp = st.plot_auc_roc_curve(
    y_true, y_hat_mass_fp, "Normalized min of atomic weight for 9-mers ROC"
)


aggrigate_atomic_weights = all_statistics.group_by("peptide").agg(
    (pl.col("9mer_weight").list.max()).first().alias("score"),
    pl.col("epitope").first().alias("epitope"),
)
y_hat_mass_fp = st.normalized_pLDDT_30mer(aggrigate_atomic_weights, "score")
y_true = aggrigate_atomic_weights.select(pl.col("epitope"))
norm_max_weight_9mer_fp = st.plot_auc_roc_curve(
    y_true, y_hat_mass_fp, "Normalized max of atomic mass for 9-mers ROC"
)
# Shows the AUC scores for each focal protein based on the mean, max, and min of the atomic weights of the 9-mers
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
