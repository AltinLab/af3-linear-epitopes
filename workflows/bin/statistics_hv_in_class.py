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
    "/scratch/sromero/af3-linear-epitopes/data/hv_class/focal_protein/staged/hv_class_focal_protein.filt.clust.parquet"
).filter(pl.col("representative"))
hv_class_dat = pl.read_parquet(
    "/scratch/sromero/af3-linear-epitopes/data/hv_class/peptide/staged/hv_class_peptide.filt.parquet"
)

all_statistics_hv_class = st.statistics(
    hv_class_dat,
    "data/hv_class/peptide/inference",
)


all_statistics_hv_class = st.peptide_9mer(
    all_statistics_hv_class,
    "data/hv_class/peptide/inference",
)


all_statistics_hv_class = st.statistics_9mer(
    all_statistics_hv_class,
    "data/hv_class/peptide/inference",
)
all_statistics_hv_class = st.pae_statistics(
    all_statistics_hv_class,
    "data/hv_class/peptide/inference",
)
all_statistics_hv_class = st.pl_avg_weight(
    all_statistics_hv_class, "data/hv_class/peptide/inference"
)
# ______________________________________________________________________
# all_statistics_hv_class = st.pl_helix(
#     all_statistics_hv_class, "data/hv_class/peptide/inference"
# )
# all_statistics_hv_class = all_statistics_hv_class.with_columns(
#     (
#         (
#             pl.col("helix").list.sum().cast(pl.Float64)
#             / pl.col("helix").list.len().cast(pl.Float64).fill_null(0)
#         )
#         * 100
#     ).alias("true__helix_percentage")
# )

# all_statistics_hv_class = st.pl_beta(
#     all_statistics_hv_class, " data/hv_class/peptide/inference"
# )
# all_statistics_hv_class = all_statistics_hv_class.with_columns(
#     (
#         (
#             pl.col("beta").list.sum().cast(pl.Float64)
#             / pl.col("beta").list.len().cast(pl.Float64).fill_null(0)
#         )
#         * 100
#     ).alias("true_beta_percentage")
# )

fp_hv_class_dat = stf.pl_fp_extract(
    fp_hv_class_dat,
    "/scratch/sromero/af3-linear-epitopes/data/hv_class/focal_protein/inference",
)
all_statistics_hv_class_fp = hv_class_dat.explode(["fp_job_names", "fp_seq_idxs"]).join(
    fp_hv_class_dat, left_on="fp_job_names", right_on="job_name"
)


all_statistics_hv_class_fp = all_statistics_hv_class_fp.explode(["fp_seq_idxs"])
all_statistics_hv_class_fp = all_statistics_hv_class_fp.with_columns(
    pl.col("pLDDT")
    .list.slice(pl.col("fp_seq_idxs"), 30)
    .list.mean()
    .alias("mean_pLDDT_slice")
)

all_statistics_hv_class_fp = stf.fp_pLDDT_score_9mer(all_statistics_hv_class_fp)
all_statistics_hv_class_fp = st.pl_9mer_weight(
    all_statistics_hv_class_fp,
    "/scratch/sromero/af3-linear-epitopes/data/hv_class/peptide/inference",
)

all_statistics_hv_class.write_parquet(
    "data/hv_class/peptide/staged/01_hv_class.features.parquet"
)
all_statistics_hv_class_fp.write_parquet(
    "data/hv_class/focal_protein/staged/01_hv_class.exploded.parquet"
)
