include { splitParquet } from 'plugin/nf-parquet'
include { MSA_WORKFLOW }from './subworkflows/tgen/af3'
include { SEQ_LIST_TO_FASTA } from './modules/tgen/af3'

workflow {

    peptide_channel = Channel.fromPath("$params.data_dir/peptide/staged/*.filt.parquet").splitParquet().map{
        row -> 
            tuple(
                [
                    id : row["job_name"],
                    protein_type : "peptide",
                ],
                [row["peptide"]],
            )
    }

    peptide_fasta_channel = SEQ_LIST_TO_FASTA(peptide_channel)

    MSA_WORKFLOW(peptide_fasta_channel)
}