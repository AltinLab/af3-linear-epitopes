
process CLEAN_IEDB_BP3 {
  label 'process_local'
  conda "envs/env.yaml"

  publishDir "$params.data_dir/$params.dset_name/focal_protein/staged", mode: 'copy'

  input:
    path raw_data

  output:
  path("*iedb_bp3*.parquet")

  script:
  """
  clean_iedb_bp3.py \
    --raw_data_path ${raw_data} \
    --discard_path iedb_bp3.discard.parquet \
    -o iedb_bp3.filt.parquet
  """
}

workflow {
  
  CLEAN_IEDB_BP3(Channel.fromPath("data/iedb_bp3/raw"))

}