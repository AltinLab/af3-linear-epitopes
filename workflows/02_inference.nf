

include { SEQ_LIST_TO_FASTA } from 'modules/tgen/af3/utils'


workflow {

    non_epitope_channel = Channel.fromPath("$params.non_epitope_path/staged/*.filt.parquet").splitParquet().map{
        row -> 
            tuple(
                [
                    job_name : row["job_name"],
                    protein_type : "peptide"
                ],
                [row["peptide"]],
            )
    }
    epitope_channel Channel.fromPath("$params.epitope_path/staged/*.filt.parquet").splitParquet().map{
        row -> 
            tuple(
                [
                    job_name : row["job_name"],
                    protein_type : "peptide"
                ],
                [row["peptide"]],
            )
    }
    focal_protein_channel = Channel.fromPath("$params.focal_protein_path/staged/*.filt.parquet").splitParquet().map{
        row -> 
            tuple(
                [
                    job_name : row["job_name"],
                    protein_type : "any"
                ],
                [row["seq"]],
            )
    }

    all_seq_channel = non_epitope_channel
        .concat(epitope_channel, focal_protein_channel)

    all_fasta_channel = SEQ_LIST_TO_FASTA(all_seq_channel)

}