params.outdir = "$params.data_dir/$params.dset_name/focal_protein"

include { PARQUET_TO_FASTA } from './modules/local/utils'

process RUN_BEPIPRED {
  queue 'gpu-v100'
  cpus '8'
  clusterOptions '--nodes=1 --ntasks=1 --gres=gpu:1 --time=1-00:00:00'
  memory '64GB'
  executor "slurm"
  tag "bp3"
  conda 'envs/bp3.yaml'
  
  input:
  path(fasta)
  
  output:
  path("*.csv")
  
  script:
  """
  export TORCH_HOME=${params.torch_home}
  
  bepipred3_CLI.py \\
      -i ${fasta} \\
      -o . \\
      -pred mjv_pred \\
      -add_seq_len \\
      -esm_dir ${params.esm_dir}
  """
}

process JOIN_BEPIPRED_INFERENCE {
    label "process_local"
    conda "envs/env.yaml"
    publishDir "${params.outdir}/staged", mode: "copy"

    input:
        path(filt_dset)
        path(bepipred_csv)

    output:
        path("*.parquet")
    
    script:
    """
    #!/usr/bin/env python

    import polars as pl
    
    filt_dset = pl.read_parquet("${filt_dset}")
    bepipred_out = pl.read_csv("${bepipred_csv}").with_columns(
        pl.col("BepiPred-3.0 linear epitope score").str.strip_chars().cast(
            pl.Float64)
        )

    bp = (
        bepipred_out.group_by("Accession", maintain_order=True).agg(
            [
                pl.col("BepiPred-3.0 score").alias("bp3_score"),
                pl.col("BepiPred-3.0 linear epitope score").alias("bp3_linear_score"),
                # pl.col("Residue").str.join().alias("validate_seq"),
            ]
        )
    ).rename({"Accession": "job_name"})

    filt_dset = filt_dset.join(bp, on="job_name")
    filt_dset.write_parquet("${filt_dset.getSimpleName()}.bp3.parquet")
    """
}

workflow {
    filt_pq = Channel.fromPath("$params.data_dir/$params.dset_name/focal_protein/staged/*.filt*.parquet")
    PARQUET_TO_FASTA(filt_pq)
    RUN_BEPIPRED(PARQUET_TO_FASTA.out)
    JOIN_BEPIPRED_INFERENCE(filt_pq, RUN_BEPIPRED.out)
}