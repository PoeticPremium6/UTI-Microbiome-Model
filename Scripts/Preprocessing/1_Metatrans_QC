#!/bin/bash

# Load necessary modules
module load miniconda3/4.7.12.1
module load python/2.7.16
source activate myjupyterlabenv

# Define directories
RAW_READS_DIR="/rUTI/Raw"
FASTQC_DIR="/rUTI/New/FASTQC"
TRIMMED_READS_DIR="/rUTI/New/Trimmed"

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
        ${RAW_READS_DIR}/${sample_id}_L001_R1_001.fastq \
        ${RAW_READS_DIR}/${sample_id}_L001_R2_001.fastq \
        -o $FASTQC_DIR -t 8
    
    # Quality trimming and adapter removal
    prinseq-lite.pl -fastq ${RAW_READS_DIR}/${sample_id}_L001_R1_001.fastq \
                    -fastq2 ${RAW_READS_DIR}/${sample_id}_L001_R2_001.fastq \
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
samples=("H25361" "H25362" "H25363" "H25364" "H25365")

# Loop through and process each sample
for sample_id in "${samples[@]}"; do
    process_sample $sample_id
done

echo "Metagenomics processing completed."
