#!/usr/bin/env python
"""
Generate the base parquet file from a protein FASTA file

Since seqs must be unique, we aggregate focal_protein_id into a list
and generate a unique ID (job_name) for each protein sequence.
"""
import polars as pl
import argparse
from af3_linear_epitopes.utils import fasta_to_polars, generate_job_name


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f",
        "--fasta_path",
        type=str,
    )
    parser.add_argument(
        "-o",
        "--output_path",
        type=str,
    )
    args = parser.parse_args()

    protein_df = (
        fasta_to_polars(args.fasta_path)
        .with_columns(
            pl.col("name").str.split("=").list.get(1).alias("raw_protein_id"),
        )
        .select("raw_protein_id", "seq")
    )

    protein_df = protein_df.group_by("seq").agg(
        pl.col("raw_protein_id").alias("raw_protein_ids"),
    )

    protein_df = (
        generate_job_name(protein_df, ["seq"], name="job_name")
        .select("job_name", "seq", "raw_protein_ids")
        .sort(by="job_name")
    )

    protein_df.write_parquet(args.output_path)
