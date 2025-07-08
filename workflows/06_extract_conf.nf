params.outdir = "$params.data_dir/$params.dset_name/focal_protein"


process EXTRACT_CONF {
    label "process_local"
    conda "envs/env.yaml"
    publishDir "${params.outdir}/staged", mode: "copy"

    input:
        path(filt_dset)

    output:
        path("*.parquet")
    
    script:
    """
    extract_conf.py \\
        --input_pq ${filt_dset} \\
        --inference_path "${params.outdir}/inference" \\
        --output "${filt_dset.getSimpleName()}.conf.parquet"
    """
}

workflow {
    filt_pq = Channel.fromPath("$params.data_dir/$params.dset_name/focal_protein/staged/*.filt*.parquet")
    EXTRACT_CONF(filt_pq)
}