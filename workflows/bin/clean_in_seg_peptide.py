#!/usr/bin/env python
"""
This script will filter original IN1 dataset into non-epitopes and epitopes using the following rules:

- Epitope peptides must have at least a z-score of >= EPITOPE_Z_SCORE_THRESH across at least EPITOPE_MIN_NUM_DONORS_THRESH focal donors
- Non-epitopes must have a z-score <= NON_EPITOPE_Z_SCORE_THRESH for all donors (not just focal donors)

Additionally, each protein is annotated with its pre-C-to-S substitution sequence
for later key joining with focal proteins
"""
import polars as pl
import argparse
from pathlib import Path
from af3_linear_epitopes.utils import generate_job_name, fasta_to_polars

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

    in1_metadata = pl.read_csv(
        root_data_path / "IN1_meta_04-09-25.tsv.gz",
        separator="\t",
        # remove PV2 (HV2) from IN1 dataset
    ).filter(pl.col("Source") != "PV2")
    in1_metadata = convert_null_species_to_int(in1_metadata)

    z_score = pl.read_csv(
        root_data_path / "IM0154_IM0162_IN1_trunc86_Z-HDI95_max.tsv.gz", separator="\t"
    )

    pre_sub_peptides = fasta_to_polars(
        root_data_path / "IN1_lcRmv.fasta", desc_as_name=True
    ).rename({"name": "FullName", "seq": "pre_sub_peptide"})

    donor_colnames = sorted(z_score.select(pl.exclude("Sequence name")).columns)

    ## this keeps indices in donor_scores_list, focal_donors, and z_score list consistent
    z_score = z_score.with_columns(pl.concat_list(donor_colnames).alias("z_score_list"))

    keep_peptides = in1_metadata

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
        pl.col("PeptideSequence").alias("peptide"),
    )

    # annotate with original seq
    out_df = out_df.join(pre_sub_peptides, on="FullName")

    out_df = (
        generate_job_name(out_df, ["peptide"], name="job_name")
        .select("job_name", "peptide", "raw_peptide_id", "epitope", "pre_sub_peptide")
        .sort(by="raw_peptide_id")
    )

    out_df.write_parquet(args.output_path)
