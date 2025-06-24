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
    "data/hv/focal_protein/staged/00_focal_protein.filt.parquet"
)
peptide_test_dat = pl.read_parquet("data/hv/peptide/staged/00_hv.filt.parquet")


all_statistics = st.statistics(
    peptide_test_dat,
    "data/hv/peptide/inference",
)


all_statistics = st.peptide_9mer(
    all_statistics,
    "data/hv/peptide/inference",
)


all_statistics = st.statistics_9mer(
    all_statistics,
    "data/hv/peptide/inference",
)
all_statistics = st.pae_statistics(
    all_statistics,
    "data/hv/peptide/inference",
)
all_statistics = st.pl_avg_weight(
    all_statistics, "/scratch/sromero/af3-linear-epitopes/data/hv/peptide/inference"
)

all_statistics = st.pl_helix(
    all_statistics, "/scratch/sromero/af3-linear-epitopes/data/hv/peptide/inference"
)
all_statistics = all_statistics.with_columns(
    (
        (
            pl.col("helix").list.sum().cast(pl.Float64)
            / pl.col("helix").list.len().cast(pl.Float64).fill_null(0)
        )
        * 100
    ).alias("true__helix_percentage")
)

all_statistics = st.pl_beta(
    all_statistics, "/scratch/sromero/af3-linear-epitopes/data/hv/peptide/inference"
)
all_statistics = all_statistics.with_columns(
    (
        (
            pl.col("beta").list.sum().cast(pl.Float64)
            / pl.col("beta").list.len().cast(pl.Float64).fill_null(0)
        )
        * 100
    ).alias("true_beta_percentage")
)

fp_test_dat = stf.pl_fp_extract(
    fp_test_dat,
    "/tgen_labs/altin/alphafold3/runs/linear_peptide/data/hv/focal_protein/inference",
)
all_statistics_fp = peptide_test_dat.explode(["fp_job_names", "fp_seq_idxs"]).join(
    fp_test_dat, left_on="fp_job_names", right_on="job_name"
)
all_statistics_fp = all_statistics_fp.explode(["fp_seq_idxs"])
all_statistics_fp = all_statistics_fp.with_columns(
    pl.col("pLDDT")
    .list.slice(pl.col("fp_seq_idxs"), 30)
    .list.mean()
    .alias("mean_pLDDT_slice")
)

all_statistics.write_parquet("data/hv/peptide/staged/01_hv.features.parquet")
all_statistics_fp.write_parquet("data/hv/peptide/staged/01_hv.exploded.parquet")
