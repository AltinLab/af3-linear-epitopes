from .sasa import PatchedSASAAnalysis
import polars as pl
import matplotlib.pyplot as plt
import numpy as np
from MDAnalysis.analysis.dssp import DSSP
from mdaf3.FeatureExtraction import *
from mdaf3.AF3OutputParser import AF3Output
from pathlib import Path
from sklearn.metrics import roc_curve, auc
import matplotlib.patches as mpatches

# adds the basic statistic data(ex: mean, min, and standard deviation) to the dataset
CHUNKSIZE = 15


def polar_charge(row):
    polar_amino_acid = ["S", "T", "N", "Q", "C", "Y", "D", "E", "K", "R", "H"]
    non_polar_amino_acid = ["A", "V", "L", "I", "P", "F", "W", "M", "G"]
    seq = row["peptide"]
    polar_count = 0
    non_polar = 0
    for i in range(0, len(seq)):
        if seq[i] in polar_amino_acid:
            polar_count += 1
        else:
            non_polar += 1
    row["polar"] = polar_count
    row["non_polar"] = non_polar
    return row


def pl_polar_charge(dataset, path):
    amino_polar_charge = split_apply_combine(dataset, polar_charge, chunksize=CHUNKSIZE)
    return amino_polar_charge


# mean rsa value
def rsa_mean(dataset):
    dataset = dataset.with_columns(
        pl.col("RSA")
        .list.slice(offset=pl.col("fp_seq_idxs"), length=30)
        .list.mean()
        .alias("mean_rsa_slice")
    )
    dataset = dataset.with_columns(
        pl.col("SA")
        .list.slice(offset=pl.col("fp_seq_idxs"), length=30)
        .list.mean()
        .alias("mean_sa_slice")
    )
    return dataset


# distance from head and tail of amino acid
def distance(dataset):
    dataset = dataset.with_columns(
        (pl.col("fp_seq_idxs") + 14.5).alias("distance_from_head")
    )
    dataset = dataset.with_columns(
        (pl.col("seq").str.len_chars() - (pl.col("fp_seq_idxs") + 14.5)).alias(
            "distance_from_tail"
        )
    )
    return dataset


# finding amino acid sequence patterns
def amino_acid_freq(row, path):
    af3 = AF3Output(Path(path) / row["job_name"])
    seq = row["peptide"]
    amino_acids = {
        "A": 0,  # Alanine
        "R": 0,  # Arginine
        "N": 0,  # Asparagine
        "D": 0,  # Aspartic Acid
        "C": 0,  # Cysteine
        "Q": 0,  # Glutamine
        "E": 0,  # Glutamic Acid
        "G": 0,  # Glycine
        "H": 0,  # Histidine
        "I": 0,  # Isoleucine
        "L": 0,  # Leucine
        "K": 0,  # Lysine
        "M": 0,  # Methionine
        "F": 0,  # Phenylalanine
        "P": 0,  # Proline
        "S": 0,  # Serine
        "T": 0,  # Threonine
        "W": 0,  # Tryptophan
        "Y": 0,  # Tyrosine
        "V": 0,  # Valine
    }
    for i in range(0, len(seq)):
        amino_acids[seq[i]] += 1
    row["most_frequent_amino_acid"] = max_key = max(amino_acids, key=amino_acids.get)
    row["amino_acid_count"] = amino_acids
    return row


def pl_amino_acids(dataset, path):
    amino_acid = split_apply_combine(
        dataset, amino_acid_freq, path, chunksize=CHUNKSIZE
    )
    return amino_acid


# finding the SASA for our datasets
def sasa_fp(row, path):
    af3 = AF3Output(Path(path) / row["job_name"])
    u = af3.get_mda_universe()
    analysis = PatchedSASAAnalysis(u)
    analysis.run()
    row["RSA"] = analysis.results.relative_residue_area[0].tolist()
    row["SA"] = analysis.results.residue_area[0].tolist()

    return row


def pl_sasa_fp(dataset, path):
    area = split_apply_combine(dataset, sasa_fp, path, chunksize=CHUNKSIZE)
    return area


# finding the avg atomic weight for the 30-mers
def avg_atomic_weight_30mer(row, path):
    af3 = AF3Output(Path(path) / row["job_name"])
    u = af3.get_mda_universe()
    weight = u.atoms.total_mass()
    avg_weight = weight
    row["atomic_weight"] = avg_weight
    return row


def avg_atomic_weight_9mer(row, path):
    af3 = AF3Output(Path(path) / row["job_name"])
    u = af3.get_mda_universe()
    index = 0
    mass = []
    for j in range(index, index + 22):
        mass.append(u.residues[j : j + 9].atoms.total_mass())
    row["9mer_weight"] = mass
    return row


def pl_9mer_weight(dataset, path):
    weight_9mer = split_apply_combine(
        dataset, avg_atomic_weight_9mer, path, chunksize=CHUNKSIZE
    )
    return weight_9mer


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


def raw_loop_indices_bool(sel):
    # find helices
    # https://docs.mdanalysis.org/2.8.0/documentation_pages/analysis/dssp.html
    beta_resindices_boolmask = DSSP(sel).run().results.dssp_ndarray[0, :, 0]
    return beta_resindices_boolmask.tolist()


# Recreating an error code:______________________________________________________
def pl_structure_error(dataset, path):
    beta_dataset = split_apply_combine(dataset, beta_error, path, chunksize=CHUNKSIZE)
    beta_helix_dataset = split_apply_combine(
        beta_dataset, helix_error, path, chunksize=CHUNKSIZE
    )
    beta_helix_loop_dataset = split_apply_combine(
        beta_helix_dataset, loop_error, path, chunksize=CHUNKSIZE
    )
    beta_helix_loop_dataset = beta_helix_loop_dataset.with_columns(
        (pl.col("helix") / 30).alias("helix_percentage"),
        (pl.col("beta") / 30).alias("beta_sheet_percentage"),
        (pl.col("loop") / 30).alias("loop_percentage"),
    )
    return beta_helix_loop_dataset


def loop_error(row, path):
    af3 = AF3Output(Path(path) / row["fp_job_names"])
    u = af3.get_mda_universe()
    row["loop"] = sum(raw_loop_indices_bool(u))
    return row


# finding the beta pleats of the 30-mer
def beta_error(row, path):
    af3 = AF3Output(Path(path) / row["fp_job_names"])
    u = af3.get_mda_universe()
    row["beta"] = sum(raw_beta_indices_bool(u))
    return row


# finding the helix's of the 30mers
def helix_error(row, path):
    af3 = AF3Output(Path(path) / row["fp_job_names"])
    u = af3.get_mda_universe()
    row["helix"] = sum(raw_helix_indices_bool(u))
    return row


# ___________________________________________________________


def pl_structure(dataset, path):
    beta_dataset = split_apply_combine(dataset, beta, path, chunksize=CHUNKSIZE)
    beta_helix_dataset = split_apply_combine(
        beta_dataset, helix, path, chunksize=CHUNKSIZE
    )
    beta_helix_loop_dataset = split_apply_combine(
        beta_helix_dataset, loop, path, chunksize=CHUNKSIZE
    )
    beta_helix_loop_dataset = beta_helix_loop_dataset.with_columns(
        (pl.col("helix") / 30).alias("helix_percentage"),
        (pl.col("beta") / 30).alias("beta_sheet_percentage"),
        (pl.col("loop") / 30).alias("loop_percentage"),
    )
    return beta_helix_loop_dataset


def pl_structure_fp(dataset, path):
    dataset = split_apply_combine(dataset, structure, path, chunksize=CHUNKSIZE)
    return dataset


def structure(row, path):
    af3 = AF3Output(Path(path) / row["fp_job_names"])
    u = af3.get_mda_universe()
    index = row["fp_seq_idxs"]
    row["loop"] = raw_loop_indices_bool(u)
    row["beta"] = raw_beta_indices_bool(u)
    row["helix"] = raw_helix_indices_bool(u)
    row["loop"] = sum(row["loop"][index : index + 30])
    row["beta"] = sum(row["beta"][index : index + 30])
    row["helix"] = sum(row["helix"][index : index + 30])
    return row


def loop(row, path):
    af3 = AF3Output(Path(path) / row["job_name"])
    u = af3.get_mda_universe()
    row["loop"] = sum(raw_loop_indices_bool(u))
    return row


# finding the beta pleats of the 30-mer
def beta(row, path):
    af3 = AF3Output(Path(path) / row["job_name"])
    u = af3.get_mda_universe()
    row["beta"] = sum(raw_beta_indices_bool(u))
    return row


def pl_beta(dataset, path):
    beta_dataset = split_apply_combine(dataset, beta, path, chunksize=CHUNKSIZE)
    return beta_dataset


# finding the helix's of the 30mers
def helix(row, path):
    af3 = AF3Output(Path(path) / row["job_name"])
    u = af3.get_mda_universe()
    row["helix"] = sum(raw_helix_indices_bool(u))
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
def normalized_pLDDT_30mer(dataset, colname: str, inverse: int):
    max_pLDDT = dataset.select(pl.col(colname)).max().item()
    print("max:" + str(max_pLDDT))
    min_pLDDT = dataset.select(pl.col(colname)).min().item()
    print("min:" + str(min_pLDDT))
    if inverse == -1:
        normalized_series = (
            dataset.with_columns(
                ((1 - (pl.col(colname) - min_pLDDT) / (max_pLDDT - min_pLDDT))).alias(
                    "normalized_pLDDT"
                )
            )
            .select(pl.col("normalized_pLDDT"))
            .to_series()
            .to_numpy()
        )
    else:
        normalized_series = (
            dataset.with_columns(
                (((pl.col(colname) - min_pLDDT) / (max_pLDDT - min_pLDDT))).alias(
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


# box and whisker plot
def display_boxplot(data, title="Box and Whisker Plot", x_label="", y_label="Value"):
    """
    Displays a box and whisker plot for the given data using Matplotlib and Polars.

    Args:
        data (list, numpy.ndarray, polars.Series, polars.DataFrame, or dict/list of such):
            The data to plot.
            - If a single list, numpy.ndarray, or polars.Series: a single box plot.
            - If a polars.DataFrame:
                - If it has one numeric column, that column will be plotted.
                - If it has multiple numeric columns, each will get a box plot.
                - If it has a 'category' column and a 'value' column, it will plot
                  box plots per category.
            - If a dictionary: keys are categories, values are lists/arrays/Series of data.
            - If a list of lists/arrays/Series: each inner list/array/Series represents a category.
        title (str, optional): The title of the plot. Defaults to "Box and Whisker Plot".
        x_label (str, optional): The label(s) for the x-axis.
                                 - If a string, used as the overall x-axis label.
                                 - If a list of strings, used as tick labels for multiple categories.
                                 Defaults to an empty string.
        y_label (str, optional): The label for the y-axis. Defaults to "Value".
    """
    plot_data = []
    category_labels = []

    # --- Data Preparation Logic using Polars ---
    if isinstance(data, (list, np.ndarray)):
        # Single dataset (list or numpy array)
        plot_data.append(data)
        category_labels.append("")  # No specific category label for a single plot
    elif isinstance(data, pl.Series):
        # Single Polars Series
        plot_data.append(data.to_list())
        category_labels.append("")
    elif isinstance(data, pl.DataFrame):
        # Polars DataFrame handling
        if "category" in data.columns and "value" in data.columns:
            # Assume long format: 'category' column for grouping, 'value' for data
            grouped = data.group_by("category").agg(
                pl.col("value").list().alias("values")
            )
            for row in grouped.iter_rows(named=True):
                plot_data.append(row["values"])
                category_labels.append(str(row["category"]))
            if not x_label:
                x_label = "Category"
        else:
            # Plot each numeric column as a separate box
            for col_name in data.columns:
                if data[col_name].dtype.is_numeric():  # Check if column is numeric
                    plot_data.append(data[col_name].to_list())
                    category_labels.append(col_name)
            if not x_label:
                x_label = "Columns"  # Default label for multiple columns

    elif isinstance(data, dict):
        # Dictionary of datasets (keys are categories)
        for key, value in data.items():
            if isinstance(value, pl.Series):
                plot_data.append(value.to_list())
            elif isinstance(value, (list, np.ndarray)):
                plot_data.append(value)
            else:
                print(
                    f"Warning: Skipping unsupported data type for key '{key}' in dictionary."
                )
                continue
            category_labels.append(str(key))
        if not x_label:
            x_label = "Category"

    elif isinstance(data, list) and all(
        isinstance(d, (list, np.ndarray, pl.Series)) for d in data
    ):
        # List of datasets (each element is a category)
        for i, dataset in enumerate(data):
            if isinstance(dataset, pl.Series):
                plot_data.append(dataset.to_list())
            elif isinstance(dataset, (list, np.ndarray)):
                plot_data.append(dataset)
            else:
                continue  # Should not happen due to all() check
            cat_label = f"Category {i+1}"
            if isinstance(x_label, list) and i < len(x_label):
                cat_label = x_label[i]
            category_labels.append(cat_label)
        if not x_label:
            x_label = "Category"
    else:
        print(
            "Error: Unsupported data format. Please provide a list, numpy array, polars Series/DataFrame, or a list/dictionary of such for multiple plots."
        )
        return

    # --- Plotting with Matplotlib ---
    if not plot_data:
        print("No valid data to plot.")
        return

    plt.figure(figsize=(8, 6))

    # Handle single vs. multiple box plots
    if (
        len(plot_data) == 1 and not category_labels[0]
    ):  # Single plot, no explicit category
        plt.boxplot(plot_data[0])
        plt.tick_params(
            axis="x", which="both", bottom=False, top=False, labelbottom=False
        )  # Hide x-axis ticks/labels
    else:
        plt.boxplot(plot_data)
        if category_labels and all(category_labels):  # If we have valid category labels
            plt.xticks(
                ticks=np.arange(1, len(category_labels) + 1),
                labels=category_labels,
                rotation=45,
                ha="right",
            )
        plt.xlabel(x_label)  # Set x-axis label if provided
    plt.title(title)
    plt.ylabel(y_label)
    plt.grid(axis="y", linestyle="--", alpha=0.7)
    plt.tight_layout()
    plt.show()


# the code below are functions to create the bar graphs and AUC curves
def plot_epitope_non_epitope_stats_9mer(
    avg_true_mean_min_9mer: float,
    avg_true_std_min_9mer: float,
    avg_false_mean_min_9mer: float,
    avg_false_std_min_9mer: float,
):
    """
    Creates a grouped bar graph comparing mean, minimum, and standard deviation
    of pLDDT values for Epitopes and Non-Epitopes.

    Args:
        avg_true_mean_min_9mer (float): Average mean pLDDT for epitopes.
        avg_true_std_min_9mer (float): Average standard deviation pLDDT for epitopes.
        avg_false_mean_min_9mer (float): Average mean pLDDT for non-epitopes.
        avg_false_std_min_9mer (float): Average standard deviation pLDDT for non-epitopes.
    """
    categories = ["Epitope", "Non-Epitope"]
    # Data for each statistic type
    mean_values = [avg_true_mean_min_9mer, avg_false_mean_min_9mer]
    std_values = [avg_true_std_min_9mer, avg_false_std_min_9mer]

    # Set up bar positions
    x = np.arange(len(categories))  # the label locations
    width = 0.25  # the width of the bars

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
    ax.set_title("pLDDT Statistics 9-mer: Epitopes vs Non-Epitopes", fontsize=16)
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
                xytext=(0, 3),  # 3 points vertical offset
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=9,
            )

    autolabel_single_bar(rects1)
    autolabel_single_bar(rects3)

    plt.tight_layout()
    return fig


def plot_epitope_non_epitope_stats_30mer(
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
    x = np.arange(len(categories))  # the label locations
    width = 0.25  # the width of the bars

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
    ax.set_title("pLDDT Statistics 30-mer: Epitopes vs Non-Epitopes", fontsize=16)
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
                xytext=(0, 3),  # 3 points vertical offset
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


# dot plot against two arrays
def plot_dot_plot(
    x_values,
    y_values,
    title="Dot Plot",
    x_label="X-axis",
    y_label="Y-axis",
    marker_style="o",
    marker_color="blue",
    alpha=0.7,
    figsize=(8, 6),
):
    """
    Plots two arrays against each other as a dot plot (scatter plot).

    Args:
        x_values (list or numpy.ndarray): The values for the x-axis.
        y_values (list or numpy.ndarray): The values for the y-axis.
                                          Must have the same length as x_values.
        title (str, optional): The title of the plot. Defaults to "Dot Plot".
        x_label (str, optional): The label for the x-axis. Defaults to "X-axis".
        y_label (str, optional): The label for the y-axis. Defaults to "Y-axis".
        marker_style (str, optional): The style of the markers. E.g., 'o' for circles,
                                      'x' for 'x's, '*' for stars. Defaults to 'o'.
        marker_color (str, optional): The color of the markers. E.g., 'blue', 'red',
                                      'green', '#FF5733'. Defaults to 'blue'.
        alpha (float, optional): The transparency of the markers (0.0 to 1.0).
                                 Useful for visualizing dense data. Defaults to 0.7.
        figsize (tuple, optional): The size of the figure (width, height) in inches.
                                   Defaults to (8, 6).
    """
    if len(x_values) != len(y_values):
        print("Error: x_values and y_values must have the same length.")
        return

    plt.figure(figsize=figsize)  # Set the figure size

    # Create the scatter plot
    plt.scatter(
        x_values, y_values, marker=marker_style, color=marker_color, alpha=alpha
    )

    # Add labels and title
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)

    plt.grid(True, linestyle="--", alpha=0.6)  # Add a grid for better readability
    plt.tight_layout()  # Adjust layout to prevent labels from overlapping
    plt.show()


def plot_dictionary_bar_chart(
    data_dict: dict[str, float],
    title: str = "Comparison of Values per Category",
    x_label: str = "Category",
    y_label: str = "Value",
    sort_by_value: bool = False,  # Set to True to sort bars by their height
):
    polar_amino_acid = ["S", "T", "N", "Q", "C", "Y", "D", "E", "K", "R", "H"]
    """
    Plots a bar chart for a dictionary where keys map to single float values.

    Args:
        data_dict (dict[str, float]): A dictionary where keys are category names (str)
                                      and values are single numerical data points (float or int).
        title (str, optional): The main title for the plot.
        x_label (str, optional): The label for the x-axis (categories).
        y_label (str, optional): The label for the y-axis (values).
        sort_by_value (bool, optional): If True, bars will be sorted by their value (height).
                                        Defaults to False (sorted by key name).
        bar_color (str, optional): The color of the bars. Defaults to 'skyblue'.
    """
    if not data_dict:
        print("Error: The input dictionary is empty. No chart to plot.")
        return

    # Extract keys and values
    if sort_by_value:
        # Sort items by value (ascending)
        sorted_items = sorted(data_dict.items(), key=lambda item: item[1])
        keys = [item[0] for item in sorted_items]
        values = [item[1] for item in sorted_items]
    else:
        # Sort items by key name (alphabetical) for consistent order if not sorting by value
        sorted_items = sorted(data_dict.items())
        keys = [item[0] for item in sorted_items]
        values = [item[1] for item in sorted_items]

    colors_by_threshold = ["blue" if v in polar_amino_acid else "red" for v in keys]
    import matplotlib.patches as mpatches

    red_patch = mpatches.Patch(color="red", label="Non-Polar amino acid")
    blue_patch = mpatches.Patch(color="blue", label="Polar amino acid")

    # Create the bar chart
    plt.figure(figsize=(10, 6))  # Adjust figure size as needed
    plt.bar(keys, values, color=colors_by_threshold)
    # Add labels and title
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)

    # Rotate x-axis labels if there are many categories to prevent overlap
    if len(keys) > 5:  # Arbitrary threshold, adjust as needed
        plt.xticks(rotation=45, ha="right")

    plt.grid(axis="y", linestyle="--", alpha=0.7)  # Add horizontal grid lines
    plt.tight_layout()  # Adjust layout to prevent labels from overlapping
    plt.legend(handles=[red_patch, blue_patch])
    plt.show()


def plot_p_values_bar_chart(
    data_dict: dict[str, float],
    title: str = "Comparison of Values per Category",
    x_label: str = "Category",
    y_label: str = "Value",
    sort_by_value: bool = False,  # Set to True to sort bars by their height
):
    polar_amino_acid = ["S", "T", "N", "Q", "C", "Y", "D", "E", "K", "R", "H"]
    """
    Plots a bar chart for a dictionary where keys map to single float values.

    Args:
        data_dict (dict[str, float]): A dictionary where keys are category names (str)
                                      and values are single numerical data points (float or int).
        title (str, optional): The main title for the plot.
        x_label (str, optional): The label for the x-axis (categories).
        y_label (str, optional): The label for the y-axis (values).
        sort_by_value (bool, optional): If True, bars will be sorted by their value (height).
                                        Defaults to False (sorted by key name).
        bar_color (str, optional): The color of the bars. Defaults to 'skyblue'.
    """
    if not data_dict:
        print("Error: The input dictionary is empty. No chart to plot.")
        return

    # Extract keys and values
    if sort_by_value:
        # Sort items by value (ascending)
        sorted_items = sorted(data_dict.items(), key=lambda item: item[1])
        keys = [item[0] for item in sorted_items]
        values = [item[1] for item in sorted_items]
    else:
        # Sort items by key name (alphabetical) for consistent order if not sorting by value
        sorted_items = sorted(data_dict.items())
        keys = [item[0] for item in sorted_items]
        values = [item[1] for item in sorted_items]
    import matplotlib.patches as mpatches

    colors_by_threshold = ["blue" if v in polar_amino_acid else "red" for v in keys]

    red_patch = mpatches.Patch(color="red", label="Non-Polar amino acid")
    blue_patch = mpatches.Patch(color="blue", label="Polar amino acid")

    # Create the bar chart
    plt.figure(figsize=(10, 6))  # Adjust figure size as needed
    plt.bar(keys, values, color=colors_by_threshold)

    # Set the y-axis to a logarithmic scale
    plt.yscale("log")

    # Add a dotted black horizontal line at y = 0.05
    # The `axhline` function adds a horizontal line across the axis.
    # `linestyle=':'` creates a dotted line, and `color='black'` sets its color.
    plt.axhline(y=0.05, color="black", linestyle=":", label="Threshold at 0.05")

    # Add labels and title
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)

    # Rotate x-axis labels if there are many categories to prevent overlap
    if len(keys) > 5:  # Arbitrary threshold, adjust as needed
        plt.xticks(rotation=45, ha="right")

    plt.grid(axis="y", linestyle="--", alpha=0.7)  # Add horizontal grid lines
    plt.tight_layout()  # Adjust layout to prevent labels from overlapping
    plt.legend(
        handles=[
            red_patch,
            blue_patch,
            mpatches.Patch(color="black", linestyle=":", label="Threshold at 0.05"),
        ],
        loc="upper left",
    )
    plt.show()
