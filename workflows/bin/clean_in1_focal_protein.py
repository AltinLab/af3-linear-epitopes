#!/usr/bin/env python
"""
Generate the base parquet file from a protein FASTA file

Since seqs must be unique, we aggregate focal_protein_id into a list
and generate a unique ID (job_name) for each protein sequence.

Many of the proteins are polyproteins- we attempt to filter these out to avoid dealing with 'in-silico cleaving'
and because AF3 has a max token length of ~5120 (https://github.com/google-deepmind/alphafold3/blob/main/docs/installation.md)

Change `MAX_PROTEIN_LENGTH` to adjust hard threshold.
"""
import polars as pl
import argparse
from af3_linear_epitopes.utils import fasta_to_polars, generate_job_name

MAX_PROTEIN_LENGTH = 1500

VALID_AA = [
    "A",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "K",
    "L",
    "M",
    "N",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "V",
    "W",
    "Y",
]

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

    protein_df = fasta_to_polars(args.fasta_path, desc_as_name=True).rename(
        {"name": "raw_protein_id"}
    )

    protein_df = protein_df.group_by("seq").agg(
        pl.col("raw_protein_id").alias("raw_protein_ids"),
    )

    protein_df = protein_df.filter(pl.col("seq").str.len_chars() <= MAX_PROTEIN_LENGTH)

    protein_df = protein_df.filter(
        ~pl.col("seq").str.contains(r"[^ACDEFGHIKLMNPQRSTVWY]")
    )

    protein_df = (
        generate_job_name(protein_df, ["seq"], name="job_name")
        .select("job_name", "seq", "raw_protein_ids")
        .sort(by="job_name")
    )

    protein_df.write_parquet(args.output_path)
