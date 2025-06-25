
process RUN_BEPIPRED {
    queue 'gpu-v100'
    cpus '8'
    clusterOptions '--nodes=1 --ntasks=1 --gres=gpu:1 --time=0:30:00'
    memory '64GB'
    executor "slurm"
    tag "bp3"
    conda 'envs/bp3.yaml'

    publishDir "${params.outdir}/bp3", mode: "copy"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*")

    script:
    """
    export TORCH_HOME=${params.torch_home}

    bepipred3_CLI.py \\
        -i ${fasta} \\
        -o . \\
        -pred mjv_pred \\
        -add_seq_len
    """
}