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
import matplotlib.pyplot as plt
import numpy as np
from MDAnalysis.analysis.dssp import DSSP
from mdaf3.FeatureExtraction import *
from mdaf3.AF3OutputParser import AF3Output
from pathlib import Path
from sklearn.metrics import roc_curve, auc


CHUNKSIZE = 15
# data_dir = (
#     "/tgen_labs/altin/alphafold3/runs/linear_peptide/data/hv/focal_protein/inference"
# )
# fp_test_dat = pl.read_parquet(
#     "/scratch/sromero/af3-linear-epitopes/data/hv/focal_protein/staged/00_focal_protein.filt.parquet"
# )


def fp_extract(row, path):

    try:
        af3 = AF3Output(Path(path) / row["job_name"])
        row["pLDDT"] = (
            af3.get_mda_universe().atoms.select_atoms("name CA").tempfactors.tolist()
        )
    except FileNotFoundError:
        row["pLDDT"] = None

    return row


def pl_fp_extract(dataset, path):
    fp_dataset = split_apply_combine(dataset, fp_extract, path, chunksize=CHUNKSIZE)
    return fp_dataset


def fp_pLDDT_score_9mer(dataset, path):
    fp_9mer_mean_pLDDT = dataset.with_columns(
        pl.col("pLDDT")
        .list.slice(pl.col("fp_seq_idxs"), 30)
        .map_elements(
            lambda x: [x[i : i + 9].mean() for i in range(0, len(x) - 8)],
            return_dtype=pl.List(pl.Float64),
        )
        .alias("pLDDT_slice_9mer")
    )
    return fp_9mer_mean_pLDDT
