
process CLEAN_CALIBER {
  label 'process_local'
  conda "envs/env.yaml"

  publishDir "$params.data_dir/$params.dset_name/focal_protein/staged", mode: 'copy'

  input:
    path raw_data

  output:
  path("*caliber*.parquet")

  script:
  """
  clean_caliber.py \
    --raw_data_path ${raw_data} \
    --discard_path 00_caliber.discard.parquet \
    -o 00_caliber.filt.parquet
  """
}

workflow {
  
  CLEAN_CALIBER(Channel.fromPath("data/caliber/raw"))

}