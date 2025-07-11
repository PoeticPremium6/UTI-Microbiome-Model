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


################
#Process rDNA
#!/bin/bash

# Load necessary modules
module load miniconda3/7.12.1
module load python/2.7.16
source activate myjupyterlabenv

# Define directories
RAW_READS_DIR="/rDNA"
FASTQC_DIR="/rDNA/FASTQC"
TRIMMED_READS_DIR="/rDNA/Trimmed"

# Ensure output directories exist
mkdir -p $FASTQC_DIR $TRIMMED_READS_DIR

# Define adapters
ADAPTER_FWD="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_REV="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

# Function to process a single sample
process_sample() {
    sample_id=$1
    echo "Processing sample: $sample_id"
    
    # FASTQC on raw reads
    /Software/FastQC/fastqc \
        ${RAW_READS_DIR}/${sample_id}_R1_rDNA.fastq \
        ${RAW_READS_DIR}/${sample_id}_R2_rDNA.fastq \
        -o $FASTQC_DIR -t 8
    
    # Quality trimming and adapter removal
    prinseq-lite.pl -fastq ${RAW_READS_DIR}/${sample_id}_R1_rDNA.fastq \
                    -fastq2 ${RAW_READS_DIR}/${sample_id}_R2_rDNA.fastq \
                    -min_len 1 -ns_max_n 5 -min_qual_mean 10 -trim_qual_left 10 \
                    -trim_qual_right 10 -trim_tail_left 10 -trim_tail_right 10 \
                    -out_good ${TRIMMED_READS_DIR}/${sample_id}_good \
                    -out_bad ${TRIMMED_READS_DIR}/${sample_id}_bad

    # FASTQC on trimmed reads
    /Software/FastQC/fastqc \
        ${TRIMMED_READS_DIR}/${sample_id}_good_1.fastq \
        ${TRIMMED_READS_DIR}/${sample_id}_good_2.fastq \
        -o $FASTQC_DIR -t 8
}

# Sample IDs
samples=("A02" "B01" "C02" "D01" "D02" "E01" "E02" "F01" "F02" "G01" "H01" "H25361" "H25362" "H25363" "H25364" "H25365")

# Loop through and process each sample
for sample_id in "${samples[@]}"; do
    process_sample $sample_id
done

echo "Metagenomics processing completed."

###########################################################################
#!/bin/bash

cd /work_beegfs/sukem113/
# Setup environment
USEARCH="./Software/usearch"
RDP_DB="/Database/rdp_16s_sp.fa"
SEQ_DIR="/rDNA"
OUT_DIR="/rDNA/Clean"
MASTER_ZOTU="${OUT_DIR}/master_zOTU.fa"

# Merge FASTQ pairs and convert to FASTA
for SAMPLE in A02 B01 C02 D01 D02 E01 E02 F01 F02 G01 H01 H25361 H25362 H25363 H25364 H25365; do
    echo "Processing $SAMPLE..."
    $USEARCH --fastq_mergepairs "${SEQ_DIR}/${SAMPLE}_R1_rDNA.fastq" --reverse "${SEQ_DIR}/${SAMPLE}_R2_rDNA.fastq" \
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
./Software/usearch --usearch_global /rDNA/Clean/master_zOTU.fa \
  --db /rDNA/Clean/master_zOTU_unique_zOTUs.fasta \
  --id 0.95 --otutabout /rDNA/Clean/ASV_counts.txt \
  --strand both


# Predict taxonomy
$USEARCH -sintax "${OUT_DIR}/master_zOTU_unique_zOTUs.fasta" -db $RDP_DB -tabbedout "${OUT_DIR}/tax.txt" -strand both -sintax_cutoff 0.45

echo "Pipeline completed."
