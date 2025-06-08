include { splitParquet } from 'plugin/nf-parquet'
include { MSA_WORKFLOW }from './subworkflows/tgen/af3'
include { SEQ_LIST_TO_FASTA } from './modules/tgen/af3'

workflow {

    focal_protein_channel = Channel.fromPath("$params.data_dir/focal_protein/staged/*.filt.parquet").splitParquet().map{
        row -> 
            tuple(
                [
                    id : row["job_name"],
                    protein_type : "any",
                ],
                [row["seq"]],
            )
    }

    focal_protein_fasta_channel = SEQ_LIST_TO_FASTA(focal_protein_channel)
 
    MSA_WORKFLOW(focal_protein_fasta_channel)
}