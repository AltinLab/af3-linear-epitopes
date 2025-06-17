// default value, should be overridden (unless testing)
params.dset_name = "test"


process CLEAN_HV1 {
  queue 'compute'
  executor "slurm"
  cpus 1
  memory '3 GB'
  clusterOptions '--nodes=1 --ntasks=1 --time=01:00:00'
  conda "envs/env.yaml"

  input:
  path hv1

  output:
  path("*.parquet")

  script:
  """
  clean_hv1.py \
    -t ${hv1} \
    --output 00_hv1.cleaned.parquet
  """
}

process CLEAN_HV2 {
  queue 'compute'
  executor "slurm"
  cpus 1
  memory '3 GB'
  clusterOptions '--nodes=1 --ntasks=1 --time=01:00:00'
  conda "envs/env.yaml"

  input:
  path hv2

  output:
  path("*.parquet")


  script:
  """
  clean_hv2.py \
    -c ${hv2} \
    -o 00_hv2.cleaned.parquet
  """
}

process COMBINE_HV {
  queue 'compute'
  executor "slurm"
  cpus 1
  memory '3 GB'
  clusterOptions '--nodes=1 --ntasks=1 --time=01:00:00'
  publishDir "$params.data_dir/$params.dset_name/peptide/staged", mode: 'copy'
  conda "envs/env.yaml"

  input:
  path hv1
  path hv2

  output:
  path("*.parquet")

  script:
  """
  #!/usr/bin/env python

  import polars as pl

  hv1 = pl.read_parquet("${hv1}").with_columns(
    pl.lit(False).alias("epitope"),
  )
  hv2 = pl.read_parquet("${hv2}").with_columns(
    pl.lit(True).alias("epitope"),
  )

  combined = pl.concat([hv1, hv2], how="vertical")

  combined.write_parquet("00_hv.cleaned.parquet")
  """
}

process CLEAN_FOCAL_PROTEIN {
  queue 'compute'
  executor "slurm"
  cpus 1
  memory '3 GB'
  clusterOptions '--nodes=1 --ntasks=1 --time=01:00:00'
  publishDir "$params.data_dir/$params.dset_name/focal_protein/staged", mode: 'copy'
  conda "envs/env.yaml"

  input:
  path fp

  output:
  path("*.parquet")

  script:
  """
  clean_focal_protein.py \
    -f ${fp} \
    -o 00_focal_protein.cleaned.parquet
  """
}

workflow {

  hv1 = CLEAN_HV1(Channel.fromPath("data/hv/raw/PV1_meta_2020-11-23.tsv"))
  hv2 = CLEAN_HV2(Channel.fromPath("data/hv/raw/HV2_annot.csv"))
  COMBINE_HV(hv1, hv2)
  
  CLEAN_FOCAL_PROTEIN(Channel.fromPath("data/hv/raw/fulldesign_2019-02-27_wGBKsw.fasta"))

}