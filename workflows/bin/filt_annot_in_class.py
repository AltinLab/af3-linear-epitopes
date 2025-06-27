#!/usr/bin/env python
"""
Take the non-epitope, epitope, and focal protein datasets, and filter them down such that:

1. Non-epitope peptides have no overlap with epitope peptides (since, i.e. HV1 contains all peptides,
    reactive or not)
2. Focal proteins only contains proteins with >= 1 30-mer from both non-epitope and epitope (as a control)
3. Non-epitope only contains peptides which are present in focal proteins after step 2
4. Epitope only contains peptides which are present in focal proteins after step 2

Then, annotate non-epitope and epitope with a list column that contains the seq_ids
of the associated focal proteins and the indices into those sequences

Since this is IN1, AA subsitutions were performed. Therefore, we join on original peptides.
"""
from af3_linear_epitopes.utils import set_30mer_indices_to_true
import polars as pl
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-p",
        "--peptide_path",
        type=str,
    )
    parser.add_argument(
        "-op",
        "--output_peptide_path",
        type=str,
    )
    parser.add_argument(
        "-f",
        "--focal_protein_path",
        type=str,
    )
    parser.add_argument(
        "-of",
        "--output_focal_protein_path",
        type=str,
    )

    args = parser.parse_args()

    p_df = pl.read_parquet(args.peptide_path)
    e_df = p_df.filter(pl.col("epitope"))
    ne_df = p_df.filter(~pl.col("epitope"))

    fp_df = pl.read_parquet(args.focal_protein_path)

    ne_df = ne_df.join(e_df.select("peptide"), on="peptide", how="anti")

    match_df_ne = (
        # for those fp sequences that contain a peptide in epitope set
        fp_df.filter(
            pl.col("seq").str.contains_any(
                ne_df.select("pre_sub_peptide").to_series().implode()
            )
        )
        # record all peptides each sequence contains
        .with_columns(
            pl.col("seq")
            .str.extract_many(
                ne_df.select("pre_sub_peptide").to_series().implode(), overlapping=True
            )
            .alias("ne_peptides")
        )
        .explode("ne_peptides")
        .rename({"ne_peptides": "ne_peptide"})
        ## temporarily convert ne_peptide to a 1-length list to allow use of find_many
        .with_columns(pl.concat_list(pl.col("ne_peptide")))
        # and find their indices
        .with_columns(
            pl.col("seq")
            .str.find_many(pl.col("ne_peptide"), overlapping=True)
            .alias("ne_seq_idx")
        )
        ## now convert tmp list back to string
        .with_columns(pl.col("ne_peptide").list.first().alias("ne_peptide"))
        .group_by("job_name", "seq")
        .agg(
            pl.col("ne_peptide").alias("ne_peptides"),
            pl.col("ne_seq_idx").alias("ne_seq_idxs"),
        )
    )

    match_df_e = (
        fp_df.filter(
            pl.col("seq").str.contains_any(
                e_df.select("pre_sub_peptide").to_series().implode()
            )
        )
        .with_columns(
            pl.col("seq")
            .str.extract_many(
                e_df.select("pre_sub_peptide").to_series().implode(), overlapping=True
            )
            .alias("e_peptides")
        )
        .explode("e_peptides")
        .rename({"e_peptides": "e_peptide"})
        .with_columns(pl.concat_list(pl.col("e_peptide")))
        .with_columns(
            pl.col("seq")
            .str.find_many(pl.col("e_peptide"), overlapping=True)
            .alias("e_seq_idx")
        )
        .with_columns(pl.col("e_peptide").list.first().alias("e_peptide"))
        .group_by("job_name", "seq")
        .agg(
            pl.col("e_peptide").alias("e_peptides"),
            pl.col("e_seq_idx").alias("e_seq_idxs"),
        )
    )

    match_df_all = match_df_e.join(match_df_ne, on=["job_name", "seq"])

    match_df_all = (
        match_df_all.with_columns(
            pl.col("seq")
            .str.split("")
            # to get boolmask with equal length to seq
            .list.eval(pl.element() == "NOT_AN_AMINO_ACID")
            .alias("epitope_30mer_segment_boolmask")
        )
        .with_columns(
            pl.col("e_seq_idxs")
            # flattens the list of indices
            .list.eval(pl.element().explode().drop_nulls()).alias("flat_e_seq_idxs"),
        )
        .with_columns(
            pl.struct(
                boolmask="epitope_30mer_segment_boolmask",
                indices="flat_e_seq_idxs",
            )
            .map_elements(set_30mer_indices_to_true, return_dtype=pl.List(pl.Boolean))
            .alias("epitope_30mer_segment_boolmask")
        )
    )

    # filter focal proteins down to only the proteins that contain >1 epitope and >=1 nonepitope
    # this is a control
    fp_df = fp_df.join(
        match_df_all.select("job_name", "epitope_30mer_segment_boolmask"), on="job_name"
    ).sort(by="job_name")

    # annotate epitope and nonepitope df with the job_names and indices into those job_names
    e_annot_filt = (
        match_df_all.select("job_name", "seq", "e_peptides", "e_seq_idxs")
        .explode("e_peptides", "e_seq_idxs")
        .rename(
            {
                "e_peptides": "pre_sub_peptide",
                "job_name": "fp_job_names",
                "e_seq_idxs": "fp_seq_idxs",
            }
        )
        .group_by("pre_sub_peptide")
        .agg(pl.col("fp_job_names"), pl.col("fp_seq_idxs"))
    )

    e_df = e_df.join(e_annot_filt, on="pre_sub_peptide")

    ne_annot_filt = (
        match_df_all.select("job_name", "seq", "ne_peptides", "ne_seq_idxs")
        .explode("ne_peptides", "ne_seq_idxs")
        .rename(
            {
                "ne_peptides": "pre_sub_peptide",
                "job_name": "fp_job_names",
                "ne_seq_idxs": "fp_seq_idxs",
            }
        )
        .group_by("pre_sub_peptide")
        .agg(pl.col("fp_job_names"), pl.col("fp_seq_idxs"))
    )

    ne_df = ne_df.join(ne_annot_filt, on="pre_sub_peptide")

    p_df = pl.concat([e_df, ne_df], how="vertical").sort(by="raw_peptide_id")

    p_df.write_parquet(args.output_peptide_path)
    fp_df.write_parquet(args.output_focal_protein_path)
