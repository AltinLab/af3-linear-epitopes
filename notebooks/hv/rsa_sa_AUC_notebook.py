
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


all_statistics_fp


all_statistics_fp.with_columns(pl.col("RSA").list.mean().alias("Mean_RSA"))


y_hat_RSA_fp = st.normalized_pLDDT_30mer(all_statistics_fp, "Mean_RSA")
y_true_RSA = all_statistics_fp.select(pl.col("epitope"))
st.plot_auc_roc_curve(
    y_true_RSA, y_hat_RSA_fp, "Normalized mean RSA values for 30mer fp ROC"
)
