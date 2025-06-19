import polars as pl
import matplotlib.pyplot as plt
import numpy as np
from MDAnalysis.analysis.dssp import DSSP
from mdaf3.FeatureExtraction import *
from mdaf3.AF3OutputParser import AF3Output
from pathlib import Path
from sklearn.metrics import roc_curve, auc

# adds the basic statistic data(ex: mean, min, and standard deviation) to the dataset
CHUNKSIZE = 15


# finding the avg atomic weight for the 30-mers
def avg_atomic_weight_30mer(row, path):
    af3 = AF3Output(Path(path) / row["job_name"])
    u = af3.get_mda_universe()
    weight = u.atoms.total_mass()
    avg_weight = weight / 30
    row["avg_atomic_weight"] = avg_weight
    return row


def pl_avg_weight(dataset, path):
    avg_weight = split_apply_combine(
        dataset, avg_atomic_weight_30mer, path, chunksize=CHUNKSIZE
    )
    return avg_weight


def raw_helix_indices_bool(sel):
    # find helices
    # https://docs.mdanalysis.org/2.8.0/documentation_pages/analysis/dssp.html
    helix_resindices_boolmask = DSSP(sel).run().results.dssp_ndarray[0, :, 1]
    return helix_resindices_boolmask.tolist()


def raw_beta_indices_bool(sel):
    # find helices
    # https://docs.mdanalysis.org/2.8.0/documentation_pages/analysis/dssp.html
    beta_resindices_boolmask = DSSP(sel).run().results.dssp_ndarray[0, :, 2]
    return beta_resindices_boolmask.tolist()


# finding the beta pleats of the 30-mer
def beta(row, path):
    af3 = AF3Output(Path(path) / row["job_name"])
    u = af3.get_mda_universe()
    row["beta"] = raw_beta_indices_bool(u)
    return row


def pl_beta(dataset, path):
    beta_dataset = split_apply_combine(dataset, beta, path, chunksize=CHUNKSIZE)
    return beta_dataset


# finding the helix's of the 30mers
def helix(row, path):
    af3 = AF3Output(Path(path) / row["job_name"])
    u = af3.get_mda_universe()

    row["helix"] = raw_helix_indices_bool(u)
    return row


def pl_helix(dataset, path):
    helix_dataset = split_apply_combine(dataset, helix, path, chunksize=CHUNKSIZE)
    return helix_dataset


# creating col of statistics(mean,min,std) of every 9mer peptide
def statistics_9mer(dataset, path):
    mean_dataset_9mer = pl_mean_9mer(dataset, path)
    min_mean_dataset_9mer = pl_min_9mer(mean_dataset_9mer, path)
    return pl_std_9mer(min_mean_dataset_9mer, path)


def peptide_9mer(dataset, path):
    col_list = []
    peptides = dataset.select("peptide").to_series()
    for j in range(0, len(peptides)):
        peptide_30mer = peptides[j]
        nine_mer_seq = []
        for i in range(0, 22):
            nine_mer_seq.append(peptide_30mer[i : i + 9])
        col_list.append(nine_mer_seq)
    dataset = dataset.with_columns(pl.Series(col_list).alias("9mer_seq"))
    return dataset


def mean_func_9mer(row, path):
    af3 = AF3Output(Path(path) / row["job_name"])
    pLDDT = af3.get_mda_universe().atoms.select_atoms("name CA").tempfactors
    pLDDT_9mer = []
    for i in range(0, 22):
        pLDDT_9mer.append(pLDDT[i : i + 9].mean())

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
    for i in range(0, 22):
        pLDDT_9mer.append(pLDDT[i : i + 9].min())

    row["9mer_min_pLDDT"] = pLDDT_9mer
    return row


def pl_min_9mer(dataset, path):
    min_dataset = split_apply_combine(dataset, min_func_9mer, path, chunksize=CHUNKSIZE)
    return min_dataset


def std_func_9mer(row, path):
    af3 = AF3Output(Path(path) / row["job_name"])
    pLDDT = af3.get_mda_universe().atoms.select_atoms("name CA").tempfactors
    pLDDT_9mer = []
    for i in range(0, 22):
        pLDDT_9mer.append(pLDDT[i : i + 9].std())

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


def statistics(dataset, path):
    mean_dataset = pl_mean(dataset, path)
    min_mean_dataset = pl_min(mean_dataset, path)
    return pl_std(min_mean_dataset, path)


# normalizing the pLDDT values
def normalized_pLDDT_30mer(dataset, colname: str):
    max_pLDDT = dataset.select(pl.col(colname)).max().item()
    print("max:" + str(max_pLDDT))
    min_pLDDT = dataset.select(pl.col(colname)).min().item()
    print("min:" + str(min_pLDDT))
    mean_pLDDT = dataset.select(pl.col(colname)).to_series()
    normalized_series = (
        dataset.with_columns(
            (1 - ((pl.col(colname) - min_pLDDT) / (max_pLDDT - min_pLDDT))).alias(
                "normalized_pLDDT"
            )
        )
        .select(pl.col("normalized_pLDDT"))
        .to_series()
        .to_numpy()
    )
    return normalized_series


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


def plot_auc_roc_curve(
    y_true: np.ndarray,
    y_scores: np.ndarray,
    title: str = "Receiver Operating Characteristic (ROC) Curve",
) -> plt.Figure:
    """
    Generates and plots the Receiver Operating Characteristic (ROC) curve and calculates
    the Area Under the Curve (AUC) for binary classification predictions.

    Args:
        y_true (np.ndarray): True binary labels (0 or 1).
        y_scores (np.ndarray): Target scores, usually the probability estimates
                                of the positive class.
        title (str, optional): Title for the plot. Defaults to
                               "Receiver Operating Characteristic (ROC) Curve".

    Returns:
        matplotlib.figure.Figure: The generated matplotlib figure object containing the ROC curve.
    """

    # Ensure inputs are numpy arrays
    y_true = np.asarray(y_true)
    y_scores = np.asarray(y_scores)

    # Calculate False Positive Rate (FPR), True Positive Rate (TPR), and thresholds
    # fpr: array of shape (n_thresholds,)
    #     Increasing false positive rates such that element i is the false positive rate
    #     of predictions with score >= thresholds[i].
    # tpr: array of shape (n_thresholds,)
    #     Increasing true positive rates such that element i is the true positive rate
    #     of predictions with score >= thresholds[i].
    # thresholds: array of shape (n_thresholds,)
    #     Decreasing thresholds on the decision function used to compute fpr and tpr.
    #     thresholds[0] represents no instances being predicted as positive.
    fpr, tpr, thresholds = roc_curve(y_true, y_scores)

    # The AUC provides an aggregate measure of performance across all possible
    # classification thresholds. It ranges from 0 to 1, where 1 is perfect
    # classification and 0.5 is random.
    roc_auc = auc(fpr, tpr)

    fig, ax = plt.subplots(figsize=(8, 8))

    ax.plot(
        fpr, tpr, color="darkorange", lw=2, label=f"ROC curve (AUC = {roc_auc:.2f})"
    )

    # A random classifier would have an AUC of 0.5, indicated by this diagonal line.
    ax.plot(
        [0, 1],
        [0, 1],
        color="navy",
        lw=2,
        linestyle="--",
        label="Random Classifier (AUC = 0.5)",
    )

    ax.set_title(title, fontsize=16)
    ax.set_xlabel("False Positive Rate (FPR)", fontsize=12)
    ax.set_ylabel("True Positive Rate (TPR)", fontsize=12)

    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.0, 1.05])

    ax.legend(loc="lower right")

    ax.grid(True, linestyle="--", alpha=0.7)

    plt.tight_layout()
    return fig
