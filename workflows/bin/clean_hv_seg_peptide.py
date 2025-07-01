#!/usr/bin/env python
"""
This script will filter original HV1 dataset into non-epitopes and epitopes using the following rules:

- Epitope peptides must have at least a z-score of >= EPITOPE_Z_SCORE_THRESH across at least EPITOPE_MIN_NUM_DONORS_THRESH donors
- Non-epitopes must have a z-score <= NON_EPITOPE_Z_SCORE_THRESH for all donors

"""
import polars as pl
import argparse
from pathlib import Path
from af3_linear_epitopes.utils import generate_job_name

EPITOPE_Z_SCORE_THRESH = 10
EPITOPE_MIN_NUM_DONORS_THRESH = 2
NON_EPITOPE_Z_SCORE_THRESH = 7

NULL_SPECIES_PLACEHOLDER = 99999999


def convert_null_species_to_int(df):
    df = df.with_columns(
        pl.when(pl.col("SpeciesID").is_not_null())
        .then(pl.col("SpeciesID"))
        .otherwise(pl.lit(NULL_SPECIES_PLACEHOLDER))
    )
    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-r",
        "--raw_data_path",
        type=str,
    )
    parser.add_argument(
        "-o",
        "--output_path",
        type=str,
    )

    args = parser.parse_args()

    root_data_path = Path(args.raw_data_path)

    hv1_metadata = pl.read_csv(
        root_data_path / "PV1_meta_2020-11-23.tsv",
        separator="\t",
    ).filter(pl.col("Category") == "SetCover")

    z_score = pl.read_csv(
        root_data_path / "SHERC_combined_wSB_6-23-21_Z-HDI95.tsv", separator="\t"
    )

    donor_colnames = z_score.select(pl.exclude("Sequence name")).columns

    z_score = z_score.with_columns(pl.concat_list(donor_colnames).alias("z_score_list"))

    keep_peptides = hv1_metadata

    # add z score lists to peptides
    keep_peptides = keep_peptides.join(
        z_score, left_on="CodeName", right_on="Sequence name"
    )
    keep_peptides = keep_peptides.with_columns(
        pl.col("z_score_list").list.drop_nulls().alias("focal_z_score_list")
    )

    keep_peptides = keep_peptides.with_columns(
        pl.col("focal_z_score_list").list.filter(
            (pl.element().is_not_null()) & (pl.element() >= EPITOPE_Z_SCORE_THRESH)
        )
    )

    keep_peptides = keep_peptides.with_columns(
        pl.col("focal_z_score_list").list.len().alias("focal_z_meet_thresh"),
        pl.col("focal_z_score_list").list.mean().alias("focal_z_mean"),
        pl.col("z_score_list").list.drop_nulls().list.max().alias("all_z_max"),
    )

    non_epitopes = keep_peptides.filter(
        pl.col("all_z_max") < NON_EPITOPE_Z_SCORE_THRESH
    ).with_columns(pl.lit(False).alias("epitope"))

    epitopes = keep_peptides.filter(
        pl.col("focal_z_meet_thresh") >= EPITOPE_MIN_NUM_DONORS_THRESH
    ).with_columns(pl.lit(True).alias("epitope"))

    out_df = pl.concat([non_epitopes, epitopes]).with_columns(
        pl.col("CodeName").alias("raw_peptide_id"),
        pl.col("Peptide").alias("peptide"),
    )

    out_df = (
        generate_job_name(out_df, ["peptide"], name="job_name")
        .select("job_name", "peptide", "raw_peptide_id", "epitope")
        .sort(by="raw_peptide_id")
    )

    out_df.write_parquet(args.output_path)
