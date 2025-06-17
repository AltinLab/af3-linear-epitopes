import polars as pl
import matplotlib.pyplot as plt
import numpy as np
from mdaf3.FeatureExtraction import *
from mdaf3.AF3OutputParser import AF3Output
from pathlib import Path

# adds the basic statistic data(ex: mean, min, and standard deviation) to the dataset
CHUNKSIZE = 15


def statistics_9mer(dataset, path):
    mean_dataset_9mer = pl_mean_9mer(dataset, path)
    min_mean_dataset_9mer = pl_min_9mer(mean_dataset_9mer, path)
    return pl_std_9mer(min_mean_dataset_9mer, path)


def peptide_9mer(dataset, path):
    col_list = []
    for j in range(0, len(dataset.select("peptide"))):
        peptide_30mer = dataset.select("peptide")[j, 0]
        nine_mer_seq = []
        for i in range(0, 23):
            nine_mer_seq.append(peptide_30mer[i : i + 8])
        col_list.append(nine_mer_seq)
    dataset = dataset.with_columns(pl.Series(col_list).alias("9mer_seq"))
    return dataset


def mean_func_9mer(row, path):
    af3 = AF3Output(Path(path) / row["job_name"])
    pLDDT = af3.get_mda_universe().atoms.select_atoms("name CA").tempfactors
    pLDDT_9mer = []
    for i in range(0, 23):
        pLDDT_9mer.append(pLDDT[i : i + 8].mean())

    row["9mer_Mean_pLDDT"] = pLDDT_9mer
    return row


def pl_mean_9mer(dataset, path):
    mean_dataset = split_apply_combine(
        dataset, mean_func_9mer, path, chunksize=CHUNKSIZE
    )
    return mean_dataset


def min_func_9mer(row, path):
    af3 = AF3Output(Path(path) / row["job_name"])
    pLDDT = af3.get_mda_universe().atoms.select_atoms("name CA").tempfactors
    pLDDT_9mer = []
    for i in range(0, 23):
        pLDDT_9mer.append(pLDDT[i : i + 8].min())

    row["9mer_min_pLDDT"] = pLDDT_9mer
    return row


def pl_min_9mer(dataset, path):
    min_dataset = split_apply_combine(dataset, min_func_9mer, path, chunksize=CHUNKSIZE)
    return min_dataset


def std_func_9mer(row, path):
    af3 = AF3Output(Path(path) / row["job_name"])
    pLDDT = af3.get_mda_universe().atoms.select_atoms("name CA").tempfactors
    pLDDT_9mer = []
    for i in range(0, 23):
        pLDDT_9mer.append(pLDDT[i : i + 8].std())

    row["9mer_std_pLDDT"] = pLDDT_9mer
    return row


def pl_std_9mer(dataset, path):
    std_dataset = split_apply_combine(dataset, std_func_9mer, path, chunksize=CHUNKSIZE)
    return std_dataset


# Now working with PAE values
def pae_statistics(dataset, path):
    mean_dataset_pae = pl_pae_mean(dataset, path)
    min_mean_dataset_pae = pl_pae_min(mean_dataset_pae, path)
    return pl_pae_std(min_mean_dataset_pae, path)


def pae_mean_func(row, path):
    af3 = AF3Output(Path(path) / row["job_name"])
    u = af3.get_mda_universe()
    residx = u.residues[0:30].resindices
    row["mean_PAE_values"] = af3.get_pae_ndarr()[residx].mean()
    return row


def pl_pae_mean(dataset, path):
    mean_dataset = split_apply_combine(
        dataset, pae_mean_func, path, chunksize=CHUNKSIZE
    )
    return mean_dataset


def pae_min_func(row, path):
    af3 = AF3Output(Path(path) / row["job_name"])
    u = af3.get_mda_universe()
    residx = u.residues[0:30].resindices
    row["min_PAE_values"] = af3.get_pae_ndarr()[residx].min()
    return row


def pl_pae_min(dataset, path):
    min_dataset = split_apply_combine(dataset, pae_min_func, path, chunksize=CHUNKSIZE)
    return min_dataset


def pae_std_func(row, path):
    af3 = AF3Output(Path(path) / row["job_name"])
    u = af3.get_mda_universe()
    residx = u.residues[0:30].resindices
    row["std_PAE_values"] = af3.get_pae_ndarr()[residx].std()
    return row


def pl_pae_std(dataset, path):
    std_dataset = split_apply_combine(dataset, pae_std_func, path, chunksize=CHUNKSIZE)
    return std_dataset


"""def bar_graph(dataset, path):
    required_cols = ["job_name", "Mean_pLDDT", "Min_pLDDT", "Std_pLDDT"]
    # Check if all required columns for plotting are present in the DataFrame
    if not all(col in dataset.columns for col in required_cols):
        missing_cols = [col for col in required_cols if col not in dataset.columns]
        print(
            f"Error: DataFrame is missing required columns for plotting: {missing_cols}."
        )
        print(
            "Please ensure your `pl_mean`, `pl_min`, and `pl_std` functions have been"
        )
        print("applied to the dataset to add these columns before calling 'bar_graph'.")
        return

    job_names = dataset["job_name"].to_list()
    mean_values = dataset["Mean_pLDDT"].to_list()
    min_values = dataset["Min_pLDDT"].to_list()
    std_values = dataset["Std_pLDDT"].to_list()

    n_jobs = len(job_names)
    if n_jobs == 0:
        print(
            "No job data to plot. The DataFrame 'dataset' is empty or contains no job names."
        )
        return

    bar_width = 0.2
    # The x locations for the groups of bars (one group for each job_name)
    index = np.arange(n_jobs)

    fig, ax = plt.subplots(figsize=(max(10, n_jobs * bar_width * 5), 7))

    bar1 = ax.bar(
        index - bar_width, mean_values, bar_width, label="Mean pLDDT", color="skyblue"
    )
    bar2 = ax.bar(
        index, min_values, bar_width, label="Minimum pLDDT", color="lightcoral"
    )
    bar3 = ax.bar(
        index + bar_width,
        std_values,
        bar_width,
        label="Standard Deviation pLDDT",
        color="lightgreen",
    )

    # Helper function to add value labels on top of each bar
    def autolabel(bars):
        for bar in bars:
            height = bar.get_height()
            ax.annotate(
                f"{height:.2f}",
                xy=(bar.get_x() + bar.get_width() / 2, height),
                xytext=(0, 3),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=9,
            )

    # Apply autolabel to each set of bars
    autolabel(bar1)
    autolabel(bar2)
    autolabel(bar3)

    # Set titles and labels for the plot
    ax.set_xlabel("Job Name", fontsize=12)
    ax.set_ylabel("pLDDT Value", fontsize=12)
    ax.set_title("pLDDT Statistics per Job", fontsize=16)
    ax.set_xticks(index)

    ax.set_xticklabels(job_names, rotation=45, ha="right", fontsize=10)
    ax.legend()
    ax.grid(axis="y", linestyle="--", alpha=0.7)

    plt.tight_layout()
    plt.show()"""


def statistics(dataset, path):
    mean_dataset = pl_mean(dataset, path)
    min_mean_dataset = pl_min(mean_dataset, path)
    return pl_std(min_mean_dataset, path)


# Is the functions that appends the mean column to the dataset
def mean_func(row, path):
    af3 = AF3Output(Path(path) / row["job_name"])
    row["Mean_pLDDT"] = (
        af3.get_mda_universe().atoms.select_atoms("name CA").tempfactors.mean()
    )
    return row


# same for mean but for minimum
def min_func(row, path):
    af3 = AF3Output(Path(path) / row["job_name"])
    row["Min_pLDDT"] = (
        af3.get_mda_universe().atoms.select_atoms("name CA").tempfactors.min()
    )
    return row


# for standard deviation
def std_func(row, path):
    af3 = AF3Output(Path(path) / row["job_name"])
    row["Std_pLDDT"] = (
        af3.get_mda_universe().atoms.select_atoms("name CA").tempfactors.std()
    )
    return row


# applies the mean function to the dataset
def pl_mean(dataset, path):
    mean_dataset = split_apply_combine(dataset, mean_func, path, chunksize=CHUNKSIZE)
    return mean_dataset


# applies the min function to the dataset
def pl_min(dataset, path):
    min_dataset = split_apply_combine(dataset, min_func, path, chunksize=CHUNKSIZE)
    return min_dataset


# applies the standard deviation to the dataset
def pl_std(dataset, path):
    std_dataset = split_apply_combine(dataset, std_func, path, chunksize=CHUNKSIZE)
    return std_dataset


if __name__ == "__main__":
    peptide_test_dat = pl.read_parquet(
        "/scratch/sromero/af3-linear-epitopes/data/test/peptide/staged/00_hv.filt.parquet"
    )

    pl_mean(
        peptide_test_dat,
        "/scratch/sromero/af3-linear-epitopes/data/test/peptide/inference",
    )
