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
        "-c",
        "--csv_path",
        type=str,
    )
    parser.add_argument(
        "-o",
        "--output_path",
        type=str,
    )
    args = parser.parse_args()

    hv2_df = (
        pl.read_csv(
            args.csv_path,
        )
        .rename(
            {
                "library_member": "raw_peptide_id",
            }
        )
        .filter(pl.col("Category") == "PV1")
    ).select("raw_peptide_id", "peptide")

    hv2_df = generate_job_name(hv2_df, ["peptide"], name="job_name").select(
        "job_name",
        "peptide",
        "raw_peptide_id",
    )

    if hv2_df.select("peptide").unique().height != hv2_df.height:
        raise ValueError("Duplicate peptides found in the HV2 dataset.")

    hv2_df.write_parquet(args.output_path)
