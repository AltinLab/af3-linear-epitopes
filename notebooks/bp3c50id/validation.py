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

# %% [markdown]
# # Load Bepipred 3 training + test data extracted from PDB
#
# "test" bool column indiciates if protein is part of BP3 test dataset (or train)
#

# %%
import polars as pl

bp3 = pl.read_parquet("../../data/bp3c50id/focal_protein/staged/bp3c50id.bp3.parquet")

bp3_test = bp3.filter(pl.col("test"))

# %%
from sklearn.metrics import roc_auc_score

# convert true/false values into one long list
y_true = (
    bp3_test.select(pl.implode("epitope_boolmask").list.eval(pl.element().explode()))
    .to_series()
    .to_list()[0]
)

y_score = (
    bp3_test.select(pl.implode("bp3_score").list.eval(pl.element().explode()))
    .to_series()
    .to_list()[0]
)

# %%
roc_auc_score(y_true, y_score)
