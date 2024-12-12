#!/bin/bash 

# Define paths
TRANSCRIPTOME="/UTI89_transcriptome.gff"
ABUNDANCE="/master_gene_abundance_UTI89.tsv"
OUTPUT_DIR="/Processed_Data"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Process transcriptome
# Extract gene ID, start, end positions, and calculate gene lengths
awk '{print $4, $5, $10}' $TRANSCRIPTOME > $OUTPUT_DIR/transcriptome_clean_step1.tsv

# Calculate gene lengths (end - start) and save to a new file
awk '{print $1, $2, $2-$1}' $OUTPUT_DIR/transcriptome_clean_step1.tsv > $OUTPUT_DIR/gene_lengths.tsv

# Remove duplicate entries and save final gene lengths
awk '!_[$0]++' $OUTPUT_DIR/gene_lengths.tsv > $OUTPUT_DIR/UTI89_Gene_Lengths.tsv

# Process gene abundance
# Sort the abundance file
sort -V $ABUNDANCE > $OUTPUT_DIR/clean_master_gene_abundance.tsv

# Step to avoid duplicate headers: Remove the first line (if already a header) and add custom headers
# Remove duplicate headers and write final result
awk 'NR==1{next} {print}' $OUTPUT_DIR/clean_master_gene_abundance.tsv | \
awk 'BEGIN {print "Gene ID\tGene Name\tA01\tA02\tB01\tB02\tC01\tC02\tD01\tD02\tE01\tE02\tF01\tF02\tG01\tH01\tH25361\tH25362\tH25363\tH25364\tH25365"} {print $0}' > $OUTPUT_DIR/clean_master_gene_abundance_UTI89.tsv

# Final message
echo "Gene abundance file processed and saved to $OUTPUT_DIR/clean_master_gene_abundance_UTI89.tsv"

#############################
import os
import pandas as pd

# Define directories
input_dir = r"\New_Gene_Abundance"
output_dir = r"\Processed_Data\UTI89"
os.makedirs(output_dir, exist_ok=True)

# Define the sample IDs
sample_ids = [
    "A01", "A02", "B01", "B02", "C01", "C02", 
    "D01", "D02", "E01", "E02", "F01", "F02", 
    "G01", "H01", "H25361", "H25362", "H25363", "H25364", "H25365"
]

# Process each sample
for sample_id in sample_ids:
    try:
        # Define input and output file paths
        input_count_file = os.path.join(
            input_dir, f"gene_abundance_GCF_000013265.1_ASM1326v1_UTI89_{sample_id}.tsv"
        )
        output_file = os.path.join(output_dir, f"gene_counts_{sample_id}.csv")
        
        # Check if the input file exists
        if not os.path.exists(input_count_file):
            print(f"Input file not found: {input_count_file}")
            continue

        # Load the input file
        df = pd.read_csv(input_count_file, sep="\t")
        
        # Add a new column for the updated gene IDs
        df["Length Difference"] = df["End"] - df["Start"]
        df["Updated Gene ID"] = df["Reference"].str.replace(".", "_") + "_" + df["Length Difference"].astype(str)
        
        # Select relevant columns and rename for clarity
        final_df = df[["Updated Gene ID", "Gene Name", "Coverage", "FPKM", "TPM"]]
        final_df.rename(columns={"Updated Gene ID": "Gene ID"}, inplace=True)
        
        # Save to CSV
        final_df.to_csv(output_file, index=False)
        print(f"Processed sample {sample_id}: output saved to {output_file}")
    
    except Exception as e:
        print(f"Error processing sample {sample_id}: {e}")

