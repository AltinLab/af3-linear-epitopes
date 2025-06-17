# ---
# jupyter:
#   jupytext:
#     custom_cell_magics: kql
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.2
#   kernelspec:
#     display_name: linear-epitope
#     language: python
#     name: python3
# ---

# %% [markdown]
# ### Overall data structure:
#
# "Focal proteins" are the original proteins where the 30-mer peptides came from
# "Peptides" refer to the 30-mer peptides that either are or are not epitopes (reactive)
#
# Inference (predicted AF3 structures and confidences) is available for both peptides and focal proteins

# %%
# %load_ext autoreload
# %autoreload 2

# %%
import polars as pl

# this one is my package
from mdaf3.AF3OutputParser import AF3Output
from mdaf3.FeatureExtraction import *
import matplotlib.pyplot as plt
import numpy as np

from pathlib import Path

# import dataframes from "staged" directory
fp_test_dat = pl.read_parquet(
    "../data/hv/focal_protein/staged/00_focal_protein.filt.parquet"
)
peptide_test_dat = pl.read_parquet("../data/hv/peptide/staged/00_hv.filt.parquet")

# %% [markdown]
# ## Dataframe structure

# %% [markdown]
# ### Focal protein dataframe
#
# Cols:
# - job_name: name of the af3 job, use it to retreive af3 output object
# - seq: the full protein sequence from which the peptides (30mers) came
# - raw_protein_id: john's original identifier for the protein (stored in list because original data contains some dupes)

# %%
fp_test_dat

# %% [markdown]
# ### Peptide dataframe
#
# Cols:
# - job_name: name of the af3 job, use it to retreive af3 output object
# - peptide: the 30mer sequence
# - raw_protein_id: john's original identifier for the peptide
# - epitope: whether or not john's assay (pepseq) detected reactivity against the peptide
# - fp_job_names: list of job_names in the focal proteins dataframe that contain this row's 30-mer peptide
# - fp_seq_idxs: at each index in the list fp_job_names, this list contains a list of integer offsets where the 30-mer peptide appears

# %%
peptide_test_dat

# %%
from af3_linear_epitopes import statistics as st

all_statistics = st.statistics(
    peptide_test_dat,
    "/tgen_labs/altin/alphafold3/runs/linear_peptide/data/hv/peptide/inference",
)


"""
st.bar_graph(
    all_statistics,
    "/scratch/sromero/af3-linear-epitopes/data/test/focal_protein/inference",
)"""

# %%
print(all_statistics)

# %%
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

# %%
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

# %%
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


# %%
def plot_epitope_non_epitope_stats(
    avg_true_mean: float,
    avg_true_min: float,
    avg_true_std: float,
    avg_false_mean: float,
    avg_false_min: float,
    avg_false_std: float,
):
    """
    Creates a grouped bar graph comparing mean, minimum, and standard deviation
    of pLDDT values for Epitopes and Non-Epitopes.

    Args:
        avg_true_mean (float): Average mean pLDDT for epitopes.
        avg_true_min (float): Average minimum pLDDT for epitopes.
        avg_true_std (float): Average standard deviation pLDDT for epitopes.
        avg_false_mean (float): Average mean pLDDT for non-epitopes.
        avg_false_min (float): Average minimum pLDDT for non-epitopes.
        avg_false_std (float): Average standard deviation pLDDT for non-epitopes.
    """
    categories = ["Epitope", "Non-Epitope"]
    # Data for each statistic type
    mean_values = [avg_true_mean, avg_false_mean]
    min_values = [avg_true_min, avg_false_min]
    std_values = [avg_true_std, avg_false_std]

    # Set up bar positions
    x = np.arange(len(categories))
    width = 0.25

    fig, ax = plt.subplots(figsize=(10, 7))

    # Create bars for Mean, Min, and Std Dev for both categories
    rects1 = ax.bar(
        x - width,
        mean_values,
        width,
        label="Mean pLDDT",
        color="skyblue",
        edgecolor="grey",
    )
    rects2 = ax.bar(
        x, min_values, width, label="Min pLDDT", color="lightcoral", edgecolor="grey"
    )
    rects3 = ax.bar(
        x + width,
        std_values,
        width,
        label="Std Dev pLDDT",
        color="lightgreen",
        edgecolor="grey",
    )

    # Add labels, title, and custom x-axis tick labels
    ax.set_ylabel("pLDDT Value", fontsize=12)
    ax.set_title("pLDDT Statistics: Epitopes vs Non-Epitopes", fontsize=16)
    ax.set_xticks(x)
    ax.set_xticklabels(categories, fontsize=12)
    ax.legend()
    ax.grid(axis="y", linestyle="--", alpha=0.7)

    # Add value labels on top of the bars
    def autolabel_single_bar(rects):
        for rect in rects:
            height = rect.get_height()
            ax.annotate(
                f"{height:.2f}",
                xy=(rect.get_x() + rect.get_width() / 2, height),
                xytext=(0, 3),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=9,
            )

    autolabel_single_bar(rects1)
    autolabel_single_bar(rects2)
    autolabel_single_bar(rects3)

    plt.tight_layout()
    return fig


pLDDT_statistics_30mer = plot_epitope_non_epitope_stats(
    avg_true_mean,
    avg_true_min,
    avg_true_std,
    avg_false_mean,
    avg_false_min,
    avg_false_std,
)
pLDDT_statistics_30mer.savefig(
    "../results/figures/pLDDT_statistics_30mer_epitope_vs_non-epitopes.png"
)

# %%
all_statistics = st.peptide_9mer(
    all_statistics,
    "/tgen_labs/altin/alphafold3/runs/linear_peptide/data/hv/peptide/inference",
)

# %%
all_statistics = st.statistics_9mer(
    all_statistics,
    "/tgen_labs/altin/alphafold3/runs/linear_peptide/data/hv/peptide/inference",
)

# %%
st.pae_statistics(
    all_statistics,
    "/tgen_labs/altin/alphafold3/runs/linear_peptide/data/hv/peptide/inference",
)

# %%
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


# %%
def plot_epitope_non_epitope_stats_9mer(
    avg_true_mean_min_9mer: float,
    avg_true_min_min_9mer: float,
    avg_true_std_min_9mer: float,
    avg_false_mean_min_9mer: float,
    avg_false_min_min_9mer: float,
    avg_false_std_min_9mer: float,
):
    """
    Creates a grouped bar graph comparing mean, minimum, and standard deviation
    of pLDDT values for Epitopes and Non-Epitopes.

    Args:
        avg_true_mean_min_9mer (float): Average mean pLDDT for epitopes.
        avg_true_min_min_9mer (float): Average minimum pLDDT for epitopes.
        avg_true_std_min_9mer (float): Average standard deviation pLDDT for epitopes.
        avg_false_mean_min_9mer (float): Average mean pLDDT for non-epitopes.
        avg_false_min_min_9mer (float): Average minimum pLDDT for non-epitopes.
        avg_false_std_min_9mer (float): Average standard deviation pLDDT for non-epitopes.
    """
    categories = ["Epitope", "Non-Epitope"]
    # Data for each statistic type
    mean_values = [avg_true_mean_min_9mer, avg_false_mean_min_9mer]
    min_values = [avg_true_min_min_9mer, avg_false_min_min_9mer]
    std_values = [avg_true_std_min_9mer, avg_false_std_min_9mer]

    # Set up bar positions
    x = np.arange(len(categories))
    width = 0.25

    fig, ax = plt.subplots(figsize=(10, 7))

    # Create bars for Mean, Min, and Std Dev for both categories
    rects1 = ax.bar(
        x - width,
        mean_values,
        width,
        label="Mean pLDDT",
        color="skyblue",
        edgecolor="grey",
    )
    rects2 = ax.bar(
        x, min_values, width, label="Min pLDDT", color="lightcoral", edgecolor="grey"
    )
    rects3 = ax.bar(
        x + width,
        std_values,
        width,
        label="Std Dev pLDDT",
        color="lightgreen",
        edgecolor="grey",
    )

    # Add labels, title, and custom x-axis tick labels
    ax.set_ylabel("pLDDT Value", fontsize=12)
    ax.set_title("pLDDT Statistics: Epitopes vs Non-Epitopes", fontsize=16)
    ax.set_xticks(x)
    ax.set_xticklabels(categories, fontsize=12)
    ax.legend()
    ax.grid(axis="y", linestyle="--", alpha=0.7)

    # Add value labels on top of the bars
    def autolabel_single_bar(rects):
        for rect in rects:
            height = rect.get_height()
            ax.annotate(
                f"{height:.2f}",
                xy=(rect.get_x() + rect.get_width() / 2, height),
                xytext=(0, 3),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=9,
            )

    autolabel_single_bar(rects1)
    autolabel_single_bar(rects2)
    autolabel_single_bar(rects3)

    plt.tight_layout()
    return fig


pLDDT_avg_9mer = plot_epitope_non_epitope_stats_9mer(
    avg_true_mean_min_9mer,
    avg_true_min_min_9mer,
    avg_true_std_min_9mer,
    avg_false_mean_min_9mer,
    avg_false_min_min_9mer,
    avg_false_std_min_9mer,
)
pLDDT_avg_9mer.savefig(
    "../results/figures/pLDDT_statistics_9mer_epitope_vs_non-epitopes.png"
)

# %%

# %%
