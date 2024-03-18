#!/bin/bash

# Setup environment
USEARCH="/Software/usearch/usearch"
RDP_DB="/home/sukem113/Database/rdp_16s_sp.fa"
SEQ_DIR="/rUTI/Final/rRNA"
OUT_DIR="/rUTI/Final/Test"
MASTER_ZOTU="${OUT_DIR}/master_zOTU.fa"

# Merge FASTQ pairs and convert to FASTA
for SAMPLE in H25361 H25362 H25363 H25364 H25365; do
    echo "Processing $SAMPLE..."
    $USEARCH --fastq_mergepairs "${SEQ_DIR}/${SAMPLE}_R1_rRNA.fastq" --reverse "${SEQ_DIR}/${SAMPLE}_R2_rRNA.fastq" \
             --fastqout "${OUT_DIR}/merged_${SAMPLE}.fq" --sample $SAMPLE
done

# Concatenate all merged files
cat ${OUT_DIR}/merged_*.fq > "${OUT_DIR}/all.merged.fq"

# Convert to FASTA
awk '(NR-1) % 4 <= 1' "${OUT_DIR}/all.merged.fq" | sed 's/@/>/' > "$MASTER_ZOTU"

# Find unique sequences
$USEARCH --fastx_uniques "$MASTER_ZOTU" --fastaout "${OUT_DIR}/master_zOTU_unique.fasta" --sizeout --relabel Uniq

# Cluster OTUs and apply chimera filtering
$USEARCH -unoise3 "${OUT_DIR}/master_zOTU_unique.fasta" -zotus "${OUT_DIR}/master_zOTU_unique_zOTUs.fasta"

# Generate count table
$USEARCH --usearch_global "$MASTER_ZOTU" --db "${OUT_DIR}/master_zOTU_unique_zOTUs.fasta" --id 0.95 \
         --otutabout "${OUT_DIR}/ASV_counts.txt"

# Predict taxonomy
$USEARCH -sintax "${OUT_DIR}/master_zOTU_unique_zOTUs.fasta" -db $RDP_DB -tabbedout "${OUT_DIR}/tax.txt" -strand both -sintax_cutoff 0.45

echo "Pipeline completed."
