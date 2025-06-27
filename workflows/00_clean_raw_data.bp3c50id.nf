
process CLEAN_IEDB_BP3 {
  label 'process_local'
  conda "envs/env.yaml"

  publishDir "$params.data_dir/$params.dset_name/focal_protein/staged", mode: 'copy'

  input:
    path raw_data

  output:
  path("*bp3c50id*.parquet")

  script:
  """
  clean_bp3c50id.py \
    --raw_data_path ${raw_data} \
    --discard_path bp3c50id.discard.parquet \
    -o bp3c50id.filt.parquet
  """
}

workflow {
  
  CLEAN_IEDB_BP3(Channel.fromPath("data/bp3c50id/raw"))

}