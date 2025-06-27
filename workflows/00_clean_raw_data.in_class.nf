include { ANNOTATE_REPRESENTATIVE_PROTEIN_WORKFLOW }from './subworkflows/local/utils'

process CLEAN_IN_CLASS {
  label "process_local"
  conda "envs/env.yaml"

  publishDir "$params.data_dir/$params.dset_name/peptide/staged", mode: 'copy'

  input:
    path raw_data

  output:
  path("*in*.parquet")

  script:
  """
  clean_in_class_peptide.py \
    -r ${raw_data} \
    -o in_class_peptide.cleaned.parquet
  """
}

process CLEAN_IN_CLASS_FOCAL_PROTEIN {
  label "process_local"
  publishDir "$params.data_dir/$params.dset_name/focal_protein/staged", mode: 'copy'
  conda "envs/env.yaml"

  input:
  path fp

  output:
  path("*.parquet")

  script:
  """
  clean_in_class_focal_protein.py \
    -f ${fp} \
    -o in_class_focal_protein.cleaned.parquet
  """
}

process FILT_AND_ANNOT_IN_CLASS {
    queue 'compute'
    executor "slurm"
    cpus 8
    memory '32 GB'
    clusterOptions '--nodes=1 --ntasks=1 --time=01:00:00'
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

  CLEAN_IN_CLASS(Channel.fromPath("data/in_class/raw"))
  
  CLEAN_IN_CLASS_FOCAL_PROTEIN(Channel.fromPath("data/in_class/raw/INF1_targets.fasta"))

  FILT_AND_ANNOT_IN_CLASS(
    CLEAN_IN_CLASS.out,
    CLEAN_IN_CLASS_FOCAL_PROTEIN.out
  )

  // this will annotate the "filt" pq with the annotated one
  ANNOTATE_REPRESENTATIVE_PROTEIN_WORKFLOW(
    FILT_AND_ANNOT_IN_CLASS.out.focal_protein
  )

}