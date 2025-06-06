#!/usr/bin/env python
"""
Generate the base parquet file from HV1 tsv
"""
import polars as pl
import argparse
from af3_linear_epitopes.utils import generate_job_name

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-t",
        "--tsv_path",
        type=str,
    )
    parser.add_argument(
        "-o",
        "--output_path",
        type=str,
    )
    args = parser.parse_args()

    hv1_df = (
        pl.read_csv(
            args.tsv_path,
            separator="\t",
            has_header=True,
        )
        .rename(
            {
                "CodeName": "raw_peptide_id",
                "Peptide": "peptide",
            }
        )
        .with_columns(
            pl.when(pl.col("FullName").str.starts_with("ID="))
            .then(
                pl.col("FullName")
                .str.split_exact(by="ID=", n=1)
                .struct.rename_fields(["tmp", "focal_protein_id"])
                .struct.field("focal_protein_id")
                .str.split(" ")
                .list.get(0)
            )
            .otherwise(pl.col("FullName"))
            .alias("focal_protein_id")
        )
        # don't keep focal protein id (for now)- we will manually find it later
    ).select("raw_peptide_id", "peptide")

    hv1_df = generate_job_name(hv1_df, ["peptide"], name="job_name").select(
        "job_name",
        "peptide",
        "raw_peptide_id",
    )

    if hv1_df.select("peptide").unique().height != hv1_df.height:
        raise ValueError("Duplicate peptides found in the HV1 dataset.")

    hv1_df.write_parquet(args.output_path)
