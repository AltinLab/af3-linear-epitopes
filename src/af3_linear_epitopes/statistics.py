import polars as pl
import matplotlib.pyplot as plt
import numpy as np
from mdaf3.FeatureExtraction import *
from mdaf3.AF3OutputParser import AF3Output
from pathlib import Path

# adds the basic statistic data(ex: mean, min, and standard deviation) to the dataset
CHUNKSIZE = 15


def statistics(dataset, path):
    Mean_dataset = pl_mean(dataset, path)
    Min_mean_dataset = pl_min(Mean_dataset, path)
    return pl_std(Min_mean_dataset, path)


def bar_graph(dataset, path):
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

    # Set up the bar positions for grouped bars
    bar_width = 0.2  # Width of each individual bar within a group
    # The x locations for the groups of bars (one group for each job_name)
    index = np.arange(n_jobs)

    # Create the figure and axes for the plot
    # Adjust figure width dynamically based on the number of jobs to prevent overcrowding
    fig, ax = plt.subplots(figsize=(max(10, n_jobs * bar_width * 5), 7))

    # Plot the bars for Mean, Minimum, and Standard Deviation
    # Bars are offset to group them visually around each x-tick (job_name)
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
            # Position the text slightly above the bar for clarity
            ax.annotate(
                f"{height:.2f}",  # Format the value to two decimal places
                xy=(bar.get_x() + bar.get_width() / 2, height),
                xytext=(0, 3),  # 3 points vertical offset from the top of the bar
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=9,
            )  # Center horizontally, align to bottom

    # Apply autolabel to each set of bars
    autolabel(bar1)
    autolabel(bar2)
    autolabel(bar3)

    # Set titles and labels for the plot
    ax.set_xlabel("Job Name", fontsize=12)
    ax.set_ylabel("pLDDT Value", fontsize=12)
    ax.set_title("pLDDT Statistics per Job", fontsize=16)
    ax.set_xticks(index)  # Set x-ticks at the center of each group
    # Rotate x-tick labels for better readability if job names are long
    ax.set_xticklabels(job_names, rotation=45, ha="right", fontsize=10)
    ax.legend()  # Display the legend for identifying Mean, Min, Std Dev bars
    ax.grid(
        axis="y", linestyle="--", alpha=0.7
    )  # Add a horizontal grid for easier value comparison

    plt.tight_layout()  # Adjust layout to prevent labels/elements from overlapping
    plt.show()  # Display the generated plot


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
