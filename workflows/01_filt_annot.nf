process FILT_AND_ANNOT {
    queue 'compute'
    executor "slurm"
    cpus 8
    memory '32 GB'
    clusterOptions '--nodes=1 --ntasks=1 --time=01:00:00'

    publishDir(
        path: {"$params.non_epitope_path/staged"},
        pattern: "${non_epitope.getSimpleName()}*",
        mode: 'copy'
    )
    publishDir(
        path: {"$params.epitope_path/staged"},
        pattern: "${epitope.getSimpleName()}*",
        mode: 'copy'
    )
    publishDir(
        path: {"$params.focal_protein_path/staged"},
        pattern: "${focal_protein.getSimpleName()}*",
        mode: 'copy'
    )

    input:
    path non_epitope
    path epitope
    path focal_protein

    output:
    path("*.filt.parquet")
    path("*.filt.parquet")
    path("*.filt.parquet")

    script:
    def non_epitope_out = non_epitope.getSimpleName() + ".filt.parquet"
    def epitope_out = epitope.getSimpleName() + ".filt.parquet"
    def focal_protein_out = focal_protein.getSimpleName() + ".filt.parquet"
    """
    filt_annot.py \
        -n ${non_epitope} \
        -on ${non_epitope_out} \
        -e ${epitope} \
        -oe ${epitope_out} \
        -f ${focal_protein} \
        -of ${focal_protein_out}
    """

}

workflow {

    FILT_AND_ANNOT(
        Channel.fromPath("$params.non_epitope_path/staged/*.cleaned.parquet"),
        Channel.fromPath("$params.epitope_path/staged/*.cleaned.parquet"),
        Channel.fromPath("$params.focal_protein_path/staged/*.cleaned.parquet")
    )

}