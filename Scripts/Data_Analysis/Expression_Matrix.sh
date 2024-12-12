#!/bin/bash

# Set directories
STRINGTIE_DIR="/stringtie-2.0.3"
ALIGNMENT_DIR="/Alignment_Output"
ANNOTATION_DIR="/New_Gene_Annotations"
TRANSCRIPTOME_DIR="/New_Transcriptome"
ABUNDANCE_DIR="/New_Gene_Abundance"
GFF_PREFIX="GCF_000013265.1_ASM1326v1_UTI89"
GFF_FILE="${ANNOTATION_DIR}/${GFF_PREFIX}.gff"

# Create output directories if they don't exist
mkdir -p "$TRANSCRIPTOME_DIR" "$ABUNDANCE_DIR"

# Run StringTie for each BAM file
for BAM_FILE in ${ALIGNMENT_DIR}/${GFF_PREFIX}_*_mRNA_cdhit_sorted.bam; do
    SAMPLE_NAME=$(basename "$BAM_FILE" | sed "s/_mRNA_cdhit_sorted.bam//")
    OUTPUT_GFF="${TRANSCRIPTOME_DIR}/${SAMPLE_NAME}.gff"
    OUTPUT_TSV="${ABUNDANCE_DIR}/gene_abundance_${SAMPLE_NAME}.tsv"
    
    # Run StringTie for FPKM (default output format is FPKM)
    "${STRINGTIE_DIR}/stringtie" "$BAM_FILE" -p 8 -G "$GFF_FILE" -e -B -o "$OUTPUT_GFF" -A "$OUTPUT_TSV"
done

# Merge all GFF files into a final transcriptome
GFF_LIST=$(ls ${TRANSCRIPTOME_DIR}/${GFF_PREFIX}_*.gff | tr '\n' ' ')
"${STRINGTIE_DIR}/stringtie" --merge -G "$GFF_FILE" -o "${TRANSCRIPTOME_DIR}/UTI89_transcriptome.gff" $GFF_LIST

# Initialize the master file with the header (Gene ID, Gene Name, and sample names)
MASTER_TSV="${ABUNDANCE_DIR}/master_gene_abundance_UTI89.tsv"
echo -e "Gene ID\tGene Name" > "$MASTER_TSV"

# Loop through each gene abundance file and extract FPKM data
for FILE in ${ABUNDANCE_DIR}/gene_abundance_${GFF_PREFIX}_*.tsv; do
    SAMPLE_NAME=$(basename "$FILE" .tsv)
    
    # Extract header (Gene ID, Gene Name, and FPKM) from the current file
    if [ ! -f "$MASTER_TSV" ]; then
        # Add sample name to the header of the master file if it doesn't already exist
        echo -e "Gene ID\tGene Name\t${SAMPLE_NAME}" >> "$MASTER_TSV"
    fi

    # Extract Gene ID, Gene Name, and FPKM values from each file
    awk -F'\t' 'NR>1 {print $1"\t"$2"\t"$NF}' "$FILE" > temp_file

    # Merge the FPKM data with the master file
    while read -r gene_id gene_name fpkm_value; do
        # Check if the gene is already in the master file
        if ! grep -q "^$gene_id" "$MASTER_TSV"; then
            # If gene is not in the master file, add it with FPKM value for this sample
            echo -e "$gene_id\t$gene_name\t$fpkm_value" >> "$MASTER_TSV"
        else
            # If gene is already in the master file, append the FPKM value for this sample
            sed -i "/^$gene_id/s/$/\t$fpkm_value/" "$MASTER_TSV"
        fi
    done < temp_file

    # Clean up temporary file
    rm temp_file
done

echo "FPKM data processed and saved successfully."

echo "StringTie processing, merging of GFF files and gene abundance files completed successfully."

# Define the input and output paths
VF_DB="/Virulence_Factor/VFDB_setB_nt.fas/VFDB_setB_nt.fas"
OUTPUT_VF="/Virulence_Factor/VFDB_setB_nt.fas/UTI89_only_virulence_factors.fas"

# Extract only lines with "Escherichia coli UTI89" entries
grep "^>.*Escherichia coli UTI89" "$VF_DB" > "$OUTPUT_VF"

echo "Filtered virulence factor database for Escherichia coli UTI89 created at $OUTPUT_VF"

# Define the input and output paths
FILTERED_VF_DB="/Virulence_Factor/VFDB_setB_nt.fas/UTI89_only_virulence_factors.fas"
OUTPUT_CSV="/Virulence_Factor/VFDB_setB_nt.fas/UTI89_gene_virulence_factors.csv"

# Drop the first and last columns
awk -F ' ' '{for(i=2; i<NF; i++) printf $i " "; print ""}' "$FILTERED_VF_DB" > "$OUTPUT_CSV"

echo "Output with first and last columns dropped saved to $OUTPUT_CSV"



