#!/usr/bin/env python
"""
Files:
- IN1_deconv_240_deconvParse_nonspecificPepsRM_maxPerSample: Enriched/reactive peptides for each donor were aggregated on a per-species basis. Each enriched peptide contributes some amount (deconv score)
    towards a species for the donors it reacted to based on the number of 7-mers the peptide shared with that species.
    Max per sample means that the highest score for each species was taken for each donor across samples taken longitudinally.
    Deconv score (cell values) can be thought of as a "likelihood that a donor was infected by a particular species"

This script will filter original IN1 dataset into IN1 (non-epitopes) and IN2 (epitopes) using the following rules:

- Epitope peptides and non-epitope peptides must have an EK thresh <= PROP_EK_THRESH and not be listed in 'removed_EK_Peptides.tsv'
- Focal species are defined as species that have a nonzero deconv (IN1_deconv...) score for >= FOCAL_SPECIES_MIN_NUM_DONORS_THRESH donors OR have a total deconv score
    across donors of >= SUM_SCORE_THRESH
- Focal donors for a focal species are defined as the donors with a non-zero deconv score for that focal species
- Epitope peptides must come from focal species
- Epitope peptides must have at least a z-score of >= EPITOPE_Z_SCORE_THRESH across at least EPITOPE_MIN_NUM_DONORS_THRESH focal donors
- Non-epitopes must come from a focal species (because they'll be filtered down to that subset anyways in the next step)
- Non-eptiopes must have a z-score <= NON_EPITOPE_Z_SCORE_THRESH for all donors (not just focal donors)

Additionally, each protein is annotated with its pre-C-to-S substitution sequence
for later key joining with focal proteins
"""
import polars as pl
import polars.selectors as cs
import argparse
from pathlib import Path
from af3_linear_epitopes.utils import generate_job_name, fasta_to_polars

FOCAL_SPECIES_MIN_NUM_DONORS_THRESH = 10
PROP_EK_THRESH = 0.401
SUM_SCORE_THRESH = 2000

EPITOPE_Z_SCORE_THRESH = 20
EPITOPE_MIN_NUM_DONORS_THRESH = 4
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
        root_data_path / "IN1_meta_04-09-25.tsv",
        separator="\t",
        # remove PV2 (HV2) from IN1 dataset
    ).filter(pl.col("Source") != "PV2")
    in1_metadata = convert_null_species_to_int(in1_metadata)

    # high in glutamate and lysine
    removed_ek_pep = pl.read_csv(
        root_data_path / "removed_EK_Peptides.tsv",
        separator="\t",
        has_header=False,
        new_columns=["CodeName"],
    )

    deconv = pl.read_csv(
        root_data_path
        / "IN1_deconv_240_deconvParse_nonspecificPepsRM_maxPerSample.tsv",
        separator="\t",
    ).select(~cs.starts_with("S-"))

    deconv = convert_null_species_to_int(deconv)

    z_score = pl.read_csv(
        root_data_path / "IM0154_IM0162_IN1_trunc86_Z-HDI95_max.tsv",
        separator="\t",
    )

    pre_sub_peptides = fasta_to_polars(
        root_data_path / "IN1_lcRmv.fasta", desc_as_name=True
    ).rename({"name": "FullName", "seq": "pre_sub_peptide"})

    donor_colnames = sorted(deconv.select(pl.exclude("SpeciesID")).columns)

    ## this keeps indices in donor_scores_list, focal_donors, and z_score list consistent
    z_score = z_score.with_columns(pl.concat_list(donor_colnames).alias("z_score_list"))

    # EK FILTERING
    in1_metadata = in1_metadata.join(removed_ek_pep, on="CodeName", how="anti")

    in1_metadata = in1_metadata.with_columns(
        (
            pl.col("PeptideSequence").str.count_matches(r"[EK]")
            / pl.col("PeptideSequence").str.len_chars()
        ).alias("prop_EK")
    ).filter(pl.col("prop_EK") <= PROP_EK_THRESH)

    # Filter peptides to only focal species

    deconv = deconv.with_columns(
        # intentionally keep nulls
        pl.concat_list(donor_colnames).alias("donor_scores_list"),
    )

    focal_species = (
        deconv.filter(
            (
                pl.col("donor_scores_list").list.drop_nulls().list.len()
                >= FOCAL_SPECIES_MIN_NUM_DONORS_THRESH
            )
            | (
                pl.col("donor_scores_list").list.drop_nulls().list.sum()
                >= SUM_SCORE_THRESH
            )
        )
        .select("SpeciesID")
        .unique()
    )

    deconv = deconv.join(focal_species, on="SpeciesID", how="inner")

    ## now keep donors with non-null, non-zero scores (to get focal donors for each species)
    deconv = deconv.with_columns(
        (
            pl.col("donor_scores_list")
            .list.eval(((pl.element().is_not_null()) & (pl.element() > 0)))
            .list.eval(pl.element().arg_true())
        ).alias("focal_donors")
    ).select("SpeciesID", "focal_donors")

    # now add focal donors to peptide df
    keep_peptides = in1_metadata.select(
        "CodeName", "FullName", "PeptideSequence", "SpeciesID"
    ).join(deconv, on="SpeciesID")

    # add z score lists to peptides
    keep_peptides = keep_peptides.join(
        z_score, left_on="CodeName", right_on="Sequence name"
    )
    keep_peptides = keep_peptides.with_columns(
        pl.col("z_score_list")
        .list.gather(pl.col("focal_donors"))
        .alias("focal_z_score_list")
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

    # REMOVED: max per species
    # keep_peptides_by_spec = keep_peptides.partition_by("SpeciesID")
    # keep_from_partitions = []
    # for species_pep in keep_peptides_by_spec:
    #     species_pep_sorted = species_pep.sort(
    #         by=["z_meet_thresh", "z_mean", "CodeName"], descending=True
    #     )
    #     species_pep_sorted = species_pep_sorted[:MAX_PEPTIDE_PER_SPECIES]
    #     keep_from_partitions.append(species_pep_sorted)

    # keep_peptides = pl.concat(keep_from_partitions)

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
