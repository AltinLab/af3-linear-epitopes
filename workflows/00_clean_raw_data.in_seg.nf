include { ANNOTATE_REPRESENTATIVE_PROTEIN_WORKFLOW }from './subworkflows/local/utils'

process CLEAN_IN_SEG {
  label "process_local"
  conda "envs/env.yaml"

  input:
    path raw_data

  output:
  path("*in*.parquet")

  script:
  """
  clean_in_seg_peptide.py \
    -r ${raw_data} \
    -o in_seg_peptide.cleaned.parquet
  """
}

process CLEAN_IN_SEG_FOCAL_PROTEIN {
  label "process_local"
  conda "envs/env.yaml"

  input:
  path fp

  output:
  path("*.parquet")

  script:
  """
  # just reuse method
  clean_in_class_focal_protein.py \
    -f ${fp} \
    -o in_seg_focal_protein.cleaned.parquet
  """
}

process FILT_AND_ANNOT_IN_SEG {
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
    # reuse script
    filt_annot_in_class.py \
        -p ${peptide} \
        -op ${peptide_out} \
        -f ${focal_protein} \
        -of ${focal_protein_out}
    """

}

// process FILT_AND_ANNOT_IN1 {
//     label "process_local"
//     conda "envs/env.yaml"

//     publishDir(
//         path: {"$params.data_dir/$params.dset_name/peptide/staged"},
//         pattern: "${peptide.getSimpleName()}*",
//         mode: 'copy'
//     )
//     publishDir(
//         path: {"$params.data_dir/$params.dset_name/focal_protein/staged"},
//         pattern: "${focal_protein.getSimpleName()}*",
//         mode: 'copy'
//     )

//     input:
//     path peptide
//     path focal_protein

//     output:
//     path("*.filt.parquet")
//     path("*.filt.parquet")

//     script:
//     def peptide_out = peptide.getSimpleName() + ".filt.parquet"
//     def focal_protein_out = focal_protein.getSimpleName() + ".filt.parquet"
//     """
//     filt_annot_in1.py \
//         -p ${peptide} \
//         -op ${peptide_out} \
//         -f ${focal_protein} \
//         -of ${focal_protein_out}
//     """

// }

workflow {

  CLEAN_IN_SEG(Channel.fromPath("data/in_seg/raw"))
  
  CLEAN_IN_SEG_FOCAL_PROTEIN(Channel.fromPath("data/in_seg/raw/INF1_targets.fasta"))

  FILT_AND_ANNOT_IN_SEG(
    CLEAN_IN_SEG.out,
    CLEAN_IN_SEG_FOCAL_PROTEIN.out
  )

  // this will annotate the "filt" pq with the annotated one
  ANNOTATE_REPRESENTATIVE_PROTEIN_WORKFLOW(
    FILT_AND_ANNOT_IN_SEG.out.focal_protein
  )

}