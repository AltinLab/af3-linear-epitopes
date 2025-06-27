include { ANNOTATE_REPRESENTATIVE_PROTEIN_WORKFLOW }from './subworkflows/local/utils'

process CLEAN_HV1_CLASS {
  label "process_local"
  conda "envs/env.yaml"

  input:
  path hv1

  output:
  path("*.parquet")

  script:
  """
  clean_hv1.py \
    -t ${hv1} \
    --output hv1_class_peptide.cleaned.parquet
  """
}

process CLEAN_HV2_CLASS {
  label "process_local"
  conda "envs/env.yaml"

  input:
  path hv2

  output:
  path("*.parquet")


  script:
  """
  clean_hv2.py \
    -c ${hv2} \
    -o hv2_class_peptide.cleaned.parquet
  """
}

process COMBINE_HV_CLASS {
  label "process_local"
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

  combined.write_parquet("hv_class_peptide.cleaned.parquet")
  """
}

process CLEAN_HV_CLASS_FOCAL_PROTEIN {
  label "process_local"
  conda "envs/env.yaml"

  input:
  path fp

  output:
  path("*.parquet")

  script:
  """
  clean_hv1_focal_protein.py \
    -f ${fp} \
    -o hv_class_focal_protein.cleaned.parquet
  """
}

process FILT_AND_ANNOT_HV_CLASS {
    label "process_local"
    conda "envs/env.yaml"

    publishDir(
        path: {"$params.data_dir/$params.dset_name/peptide/staged"},
        pattern: "*peptide*",
        mode: 'copy'
    )
    // publishDir(
    //     path: {"$params.data_dir/$params.dset_name/focal_protein/staged"},
    //     pattern: "*focal_protein*",
    //     mode: 'copy'
    // )

    input:
    path peptide
    path focal_protein

    output:
    path("*focal_protein*.parquet"), emit: focal_protein
    path("*peptide*.parquet"), emit: peptide

    script:
    def peptide_out = peptide.getSimpleName() + ".filt.parquet"
    def focal_protein_out = focal_protein.getSimpleName() + ".filt.parquet"
    """
    filt_annot_hv_class.py \
        -p ${peptide} \
        -op ${peptide_out} \
        -f ${focal_protein} \
        -of ${focal_protein_out}
    """

}

workflow {

  hv1 = CLEAN_HV1_CLASS(Channel.fromPath("data/hv_class/raw/PV1_meta_2020-11-23.tsv"))
  hv2 = CLEAN_HV2_CLASS(Channel.fromPath("data/hv_class/raw/HV2_annot.csv"))
  COMBINE_HV_CLASS(hv1, hv2)
  
  CLEAN_HV_CLASS_FOCAL_PROTEIN(Channel.fromPath("data/hv_class/raw/fulldesign_2019-02-27_wGBKsw.fasta"))

  FILT_AND_ANNOT_HV_CLASS(COMBINE_HV_CLASS.out, CLEAN_HV_CLASS_FOCAL_PROTEIN.out)

  // this will annotate the "filt" pq with the annotated one
  ANNOTATE_REPRESENTATIVE_PROTEIN_WORKFLOW(
    FILT_AND_ANNOT_HV_CLASS.out.focal_protein
  )
}