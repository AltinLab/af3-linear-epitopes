#!/usr/bin/env python
from mdaf3.FeatureExtraction import split_apply_combine
from af3_linear_epitopes.af3_feat import extract_residue_pLDDT
from pathlib import Path
import polars as pl
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input_pq",
        type=str,
    )
    parser.add_argument(
        "--inference_path",
        type=str,
    )
    parser.add_argument(
        "--output",
        type=str,
    )

    args = parser.parse_args()

    fp = pl.read_parquet(args.input_pq)

    fp = split_apply_combine(
        fp, extract_residue_pLDDT, Path(args.inference_path), chunksize=15
    )

    fp.write_parquet(args.output)
