#!/bin/bash

# Paths
SORTMERNA_EXEC="/Software/samsa2/programs/sortmerna/sortmerna"
INDEXDB_RNA_EXEC="/Software/samsa2/programs/sortmerna/indexdb_rna"
MERGE_SCRIPT="/Software/samsa2/programs/sortmerna/scripts/merge-paired-reads.sh"
UNMERGE_SCRIPT="/Software/samsa2/programs/sortmerna/scripts/unmerge-paired-reads.sh"
RAW_READS_DIR="/rUTI/FASTA"
INDEX_DIR="/Software/samsa2/programs/sortmerna/index"
DB_DIR="/Software/samsa2/programs/sortmerna/rRNA_databases"
ALIGNED_READS_DIR="/rUTI/Aligned_Reads"
SORTMERNA_OUT_DIR="/rUTI/SortmeRNA"

# Ensure output directories exist
mkdir -p $INDEX_DIR $ALIGNED_READS_DIR $SORTMERNA_OUT_DIR

# Indexing rRNA databases
echo "Indexing rRNA databases..."
sortmerna_index() {
    local ref=$1
    local index_name=$2
    $INDEXDB_RNA_EXEC --ref ${DB_DIR}/${ref},${INDEX_DIR}/${index_name} -v
}

# List of databases to index
sortmerna_index "silva-bac-16s-id90.fasta" "silva-bac-16s-db"
sortmerna_index "silva-bac-23s-id98.fasta" "silva-bac-23s-db"
# Add other databases as needed

# Process samples
echo "Processing samples..."
process_sample() {
    local sample_id=$1
    echo "Processing ${sample_id}..."
    
    # Merge paired reads
    $MERGE_SCRIPT ${RAW_READS_DIR}/${sample_id}_R1.fa ${RAW_READS_DIR}/${sample_id}_R2.fa ${SORTMERNA_OUT_DIR}/${sample_id}_merged.fastq

    # Filter rRNA
    $SORTMERNA_EXEC --fastx -a 16 --paired_out \
    --ref ${DB_DIR}/silva-bac-16s-id90.fasta,${INDEX_DIR}/silva-bac-16s-db:\
          ${DB_DIR}/silva-bac-23s-id98.fasta,${INDEX_DIR}/silva-bac-23s-db \
    --reads ${SORTMERNA_OUT_DIR}/${sample_id}_merged.fastq \
    --other ${ALIGNED_READS_DIR}/${sample_id}_other \
    --aligned ${ALIGNED_READS_DIR}/${sample_id}_aligned \
    --log -v

    # Unmerge rRNA and mRNA reads
    $UNMERGE_SCRIPT ${ALIGNED_READS_DIR}/${sample_id}_aligned.fastq ${ALIGNED_READS_DIR}/${sample_id}_R1_rRNA.fastq ${ALIGNED_READS_DIR}/${sample_id}_R2_rRNA.fastq
    $UNMERGE_SCRIPT ${ALIGNED_READS_DIR}/${sample_id}_other.fastq ${ALIGNED_READS_DIR}/${sample_id}_R1_mRNA.fastq ${ALIGNED_READS_DIR}/${sample_id}_R2_mRNA.fastq
}

# Sample IDs to process
samples=("H25361" "H25362" "H25363" "H25364" "H25365")

# Loop through and process each sample
for sample_id in "${samples[@]}"; do
    process_sample $sample_id
done

echo "rRNA filtering completed."
