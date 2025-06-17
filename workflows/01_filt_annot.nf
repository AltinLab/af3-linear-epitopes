process FILT_AND_ANNOT {
    queue 'compute'
    executor "slurm"
    cpus 8
    memory '32 GB'
    clusterOptions '--nodes=1 --ntasks=1 --time=01:00:00'
    conda "envs/env.yaml"

    publishDir(
        path: {"$params.data_dir/$params.dset_name/peptide/staged"},
        pattern: "${peptide.getSimpleName()}*",
        mode: 'copy'
    )
    publishDir(
        path: {"$params.data_dir/$params.dset_name/focal_protein/staged"},
        pattern: "${focal_protein.getSimpleName()}*",
        mode: 'copy'
    )

    input:
    path peptide
    path focal_protein

    output:
    path("*.filt.parquet")
    path("*.filt.parquet")

    script:
    def peptide_out = peptide.getSimpleName() + ".filt.parquet"
    def focal_protein_out = focal_protein.getSimpleName() + ".filt.parquet"
    """
    filt_annot.py \
        -p ${peptide} \
        -op ${peptide_out} \
        -f ${focal_protein} \
        -of ${focal_protein_out}
    """

}

workflow {

    FILT_AND_ANNOT(
        Channel.fromPath("$params.data_dir/$params.dset_name/peptide/staged/*.cleaned.parquet"),
        Channel.fromPath("$params.data_dir/$params.dset_name/focal_protein/staged/*.cleaned.parquet"),
    )

}