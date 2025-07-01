include { ANNOTATE_REPRESENTATIVE_PROTEIN_WORKFLOW }from './subworkflows/local/utils'

process CLEAN_HV_SEG {
  label "process_local"
  conda "envs/env.yaml"

  input:
    path raw_data

  output:
  path("*hv_seg*.parquet")

  script:
  """
  clean_hv_seg_peptide.py \\
    -r ${raw_data} \\
    -o hv_seg_peptide.cleaned.parquet
  """
}

process CLEAN_HV_SEG_FOCAL_PROTEIN {
  label "process_local"
  conda "envs/env.yaml"

  input:
  path fp

  output:
  path("*.parquet")

  script:
  """
  clean_hv1_focal_protein.py \\
    -f ${fp} \\
    -o hv_seg_focal_protein.cleaned.parquet
  """
}

process FILT_AND_ANNOT_HV_SEG {
    label "process_local"
    conda "envs/env.yaml"

    // publishDir(
    //     path: {"$params.data_dir/$params.dset_name/peptide/staged"},
    //     pattern: "*peptide*",
    //     mode: 'copy'
    // )
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
    # just reuse this script
    filt_annot_hv_class.py \\
        -p ${peptide} \\
        -op ${peptide_out} \\
        -f ${focal_protein} \\
        -of ${focal_protein_out}
    """

}

workflow {

  CLEAN_HV_SEG(Channel.fromPath("data/hv_seg/raw"))
  CLEAN_HV_SEG_FOCAL_PROTEIN(Channel.fromPath("data/hv_seg/raw/fulldesign_2019-02-27_wGBKsw.fasta"))

  FILT_AND_ANNOT_HV_SEG(CLEAN_HV_SEG.out, CLEAN_HV_SEG_FOCAL_PROTEIN.out)

  // this will annotate the "filt" pq with the annotated one
  ANNOTATE_REPRESENTATIVE_PROTEIN_WORKFLOW(
    FILT_AND_ANNOT_HV_SEG.out.focal_protein
  )
}