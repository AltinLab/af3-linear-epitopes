from pathlib import Path
import polars as pl


FORMAT_30MER_COLS = [
    "code_name",
    "peptide",
    "category",
    "species_id",
]

# class MyBigDataset:
#     def __init__(self, folder: str | Path):
#         self.folder = Path(folder)
#         # gather all parquet paths
#         self._files = list(self.folder.glob("*.parquet"))
#         self._df_cache = None

#     def _load_all(self) -> pl.DataFrame:
#         # e.g. use scan_parquet for efficiency:
#         return pl.scan_parquet([str(p) for p in self._files]).collect()

#     @property
#     def df(self) -> pl.DataFrame:
#         if self._df_cache is None:
#             self._df_cache = self._load_all()
#         return self._df_cache

#     def filter_by_mhc_class(self, mhc_class: str) -> pl.DataFrame:
#         return self.df.filter(pl.col("mhc_class") == mhc_class)

#     def group_summary(self, group_cols: list[str]) -> pl.DataFrame:
#         return self.df.groupby(group_cols).agg(pl.all().mean().suffix("_mean"))
