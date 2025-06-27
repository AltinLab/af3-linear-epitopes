include { PARQUET_TO_FASTA;
            FASTA_TO_PARQUET;
            CLUSTER_FASTA;
            ANNOTATE_REPRESENTATIVES;} from '../../../modules/local/utils'

workflow ANNOTATE_REPRESENTATIVE_PROTEIN_WORKFLOW {
    take:
    focal_protein_parquet

    main:
    PARQUET_TO_FASTA(focal_protein_parquet)
    CLUSTER_FASTA(PARQUET_TO_FASTA.out)
    FASTA_TO_PARQUET(CLUSTER_FASTA.out)
    ANNOTATE_REPRESENTATIVES(focal_protein_parquet, FASTA_TO_PARQUET.out)
}