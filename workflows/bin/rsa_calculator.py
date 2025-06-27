#!/usr/bin/env python
import polars as pl
from af3_linear_epitopes import statistics as st
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-pq",
        "--input_parquet",
        type=str,
    )
    parser.add_argument(
        "-i",
        "--inference_path",
        type=str,
    )
    parser.add_argument(
        "-o",
        "--output_path",
        type=str,
    )
    args = parser.parse_args()

    all_statistics_fp = pl.read_parquet(args.input_parquet)
    all_statistics_fp = st.pl_sasa_fp(all_statistics_fp, args.inference_path)

    all_statistics_fp.write_parquet(args.output_path)
