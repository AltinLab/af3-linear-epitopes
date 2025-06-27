process PARQUET_TO_FASTA {
  label "process_local"
  conda "envs/env.yaml"
  
  input:
      path(parquet)
  
  output:
      path("*.fasta")
  
  script:
  """
  #!/usr/bin/env python

  import polars as pl
  
  df = pl.read_parquet("${parquet}")
  
  with open("from_pq.fasta", "w") as f:
      for row in df.iter_rows(named=True):
          f.write(f">{row['job_name']}\\n{row['seq']}\\n")
  """
}

process FASTA_TO_PARQUET {
  label "process_local"
  conda "envs/env.yaml"
  
  input:
      path(fasta)
  
  output:
      path("*.parquet")
  
  script:
  """
  #!/usr/bin/env python

  import polars as pl
  from af3_linear_epitopes.utils import fasta_to_polars
  
  df = fasta_to_polars("${fasta}").rename({"name" : "job_name"})
  
  df.write_parquet("${fasta.getSimpleName()}.parquet")
  """
}


process CLUSTER_FASTA {
    /*

    */
    label "process_local"
    conda "envs/mmseqs2.yaml"

    input:
        path(fasta)

    output:
        path("*.fasta")

    script:
    """
    mmseqs createdb ${fasta} DB && \\
    mmseqs cluster DB DB_clu tmp \\
        --min-seq-id 0.7 && \\
    mmseqs createsubdb DB_clu DB DB_clu_rep && \\
    mmseqs convert2fasta DB_clu_rep "${fasta.getSimpleName()}.clust.fasta"
    """
}


process ANNOTATE_REPRESENTATIVES {
    label "process_local"
    conda "envs/env.yaml"
    publishDir "$params.data_dir/$params.dset_name/focal_protein/staged", mode: 'copy'

    input:
        path(orig_pq)
        path(rep_pq)

    output:
        path("*.parquet")

    script:
    """
    #!/usr/bin/env python

    import polars as pl

    orig_pq = pl.read_parquet("${orig_pq}")

    rep_pq = pl.read_parquet("${rep_pq}").with_columns(pl.lit(True).alias("representative"))

    out_df = orig_pq.join(rep_pq.select(pl.exclude("seq")), on="job_name", how="left").with_columns(
        pl.when(pl.col("representative").is_null())
        .then(pl.lit(False))
        .otherwise(pl.col("representative"))
        .alias("representative")
    )

    out_df.write_parquet("${orig_pq.getSimpleName()}.filt.clust.parquet")
    """
}