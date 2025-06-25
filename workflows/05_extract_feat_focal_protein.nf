process EXTRACT_RSA {
    queue 'compute'
    cpus '8'
    clusterOptions '--time=1-00:00:00'
    memory '64GB'
    executor "slurm"
    tag "rsa"
    conda 'envs/env.yaml'
    publishDir "$params.data_dir/$params.dset_name/focal_protein/staged", mode: 'copy'
  
    input:
    path(pq)
    path(inf)

    output:
    path("*.parquet")

    script:
    """
    rsa_calculator.py \\
        -pq ${pq} \\
        -i ${inf} \\
        -o "05_focal_protein.rsa.parquet"
    """
}

workflow {
    filt_pq = Channel.fromPath("$params.data_dir/$params.dset_name/focal_protein/staged/*.filt.parquet")
    inf = Channel.fromPath("$params.data_dir/$params.dset_name/focal_protein/inference")

    EXTRACT_RSA(
        filt_pq,
        inf
    )
}