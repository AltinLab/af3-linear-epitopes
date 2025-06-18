import polars as pl

# this one is my package
from mdaf3.AF3OutputParser import AF3Output
from mdaf3.FeatureExtraction import *
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import roc_curve, auc
from pathlib import Path
from af3_linear_epitopes import statistics as st

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

all_statistics.write_parquet("data/hv/peptide/staged/01_hv.features.parquet")
