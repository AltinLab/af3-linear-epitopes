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
# all_statistics = pl.read_parquet(
#     "/scratch/sromero/af3-linear-epitopes/data/hv/peptide/staged/01_hv.features.parquet"
# )
all_statistics_fp = pl.read_parquet(
    "/scratch/sromero/af3-linear-epitopes/data/hv/peptide/staged/01_hv.exploded.parquet"
)
all_statistics_fp = st.pl_sasa_fp(
    all_statistics_fp,
    "/scratch/sromero/af3-linear-epitopes/data/hv/peptide/inference",
)


# all_statistics.write_parquet("data/hv/peptide/staged/01_hv.features.parquet")
all_statistics_fp.write_parquet("data/hv/peptide/staged/01_hv.exploded.parquet")
