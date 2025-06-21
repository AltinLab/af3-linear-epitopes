#!/usr/bin/env python
"""
Generate the base parquet file from CALIBER paper data
"""
import polars as pl
from pathlib import Path
import argparse
from af3_linear_epitopes.utils import generate_job_name, fasta_to_polars

MAX_PROTEIN_LENGTH = 1500


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--raw_data_path",
        type=str,
    )
    parser.add_argument(
        "--discard_path",
        type=str,
    )
    parser.add_argument(
        "-o",
        "--output_path",
        type=str,
    )
    args = parser.parse_args()

    root_data_path = Path(args.raw_data_path)

    test_df = (
        fasta_to_polars(root_data_path / "linear_test.csv")
        .rename({"name": "raw_protein_id"})
        .with_columns(pl.lit(True).alias("test"))
    )
    train_df = (
        fasta_to_polars(root_data_path / "linear_train.csv")
        .rename({"name": "raw_protein_id"})
        .with_columns(pl.lit(False).alias("test"))
    )

    caliber_df = pl.concat([test_df, train_df])

    caliber_df = caliber_df.with_columns(
        pl.col("seq")
        .str.split("")
        .list.eval(pl.element().str.contains(r"[ACDEFGHIKLMNPQRSTVWY]"))
        .alias("epitope_boolmask")
    ).with_columns(pl.col("seq").str.to_uppercase().alias("seq"))

    bad_seq = caliber_df.filter(pl.col("seq").str.contains(r"[^ACDEFGHIKLMNPQRSTVWY]"))

    if bad_seq.height > 0:
        raise ValueError(
            f"Found {bad_seq.height} sequences with invalid characters\n {bad_seq}"
        )

    caliber_df = generate_job_name(
        caliber_df,
        ["seq"],
        name="job_name",
    ).select(
        "job_name",
        "seq",
        "epitope_boolmask",
        "test",
        "raw_protein_id",
    )

    discard = caliber_df.filter(pl.col("seq").str.len_chars() > MAX_PROTEIN_LENGTH)

    caliber_df = caliber_df.filter(pl.col("seq").str.len_chars() <= MAX_PROTEIN_LENGTH)

    discard.write_parquet(
        args.discard_path,
    )

    caliber_df.write_parquet(
        args.output_path,
    )
