#!/bin/bash

# Set directories
STRINGTIE_DIR="/work_beegfs/sukem113/Software/stringtie-2.0.3"
ALIGNMENT_DIR="/work_beegfs/sukem113/rUTI/Final/Alignment_Output"
ANNOTATION_DIR="/work_beegfs/sukem113/rUTI/Final/New_Gene_Annotations"
TRANSCRIPTOME_DIR="/work_beegfs/sukem113/rUTI/Final/New_Transcriptome"
ABUNDANCE_DIR="/work_beegfs/sukem113/rUTI/Final/New_Gene_Abundance"
GFF_PREFIX="GCF_000013265.1_ASM1326v1_UTI89"
GFF_FILE="${ANNOTATION_DIR}/${GFF_PREFIX}.gff"

# Create output directories if they don't exist
mkdir -p "$TRANSCRIPTOME_DIR" "$ABUNDANCE_DIR"

# Run StringTie for each BAM file
for BAM_FILE in ${ALIGNMENT_DIR}/${GFF_PREFIX}_*_mRNA_cdhit_sorted.bam; do
    SAMPLE_NAME=$(basename "$BAM_FILE" | sed "s/_mRNA_cdhit_sorted.bam//")
    OUTPUT_GFF="${TRANSCRIPTOME_DIR}/${SAMPLE_NAME}.gff"
    OUTPUT_TSV="${ABUNDANCE_DIR}/gene_abundance_${SAMPLE_NAME}.tsv"
    
    "${STRINGTIE_DIR}/stringtie" "$BAM_FILE" -p 8 -G "$GFF_FILE" -e -B -o "$OUTPUT_GFF" -A "$OUTPUT_TSV"
done

# Merge all GFF files into a final transcriptome
GFF_LIST=$(ls ${TRANSCRIPTOME_DIR}/${GFF_PREFIX}_*.gff | tr '\n' ' ')
"${STRINGTIE_DIR}/stringtie" --merge -G "$GFF_FILE" -o "${TRANSCRIPTOME_DIR}/UTI89_transcriptome.gff" $GFF_LIST

# Initialize the master file with the header (Gene ID, Gene Name, and sample names)
MASTER_TSV="${ABUNDANCE_DIR}/master_gene_abundance_UTI89.tsv"
echo -e "Gene ID\tGene Name" > "$MASTER_TSV"

# Loop through each gene abundance file and extract TPM data
for FILE in ${ABUNDANCE_DIR}/gene_abundance_${GFF_PREFIX}_*.tsv; do
    SAMPLE_NAME=$(basename "$FILE" .tsv)
    
    # Extract header (Gene ID, Gene Name, and TPM) from the current file
    if [ ! -f "$MASTER_TSV" ]; then
        # Add sample name to the header of the master file if it doesn't already exist
        echo -e "Gene ID\tGene Name\t${SAMPLE_NAME}" >> "$MASTER_TSV"
    fi

    # Extract Gene ID, Gene Name, and TPM values from each file
    awk -F'\t' 'NR>1 {print $1"\t"$2"\t"$NF}' "$FILE" > temp_file

    # Merge the TPM data with the master file
    while read -r gene_id gene_name tpm_value; do
        # Check if the gene is already in the master file
        if ! grep -q "^$gene_id" "$MASTER_TSV"; then
            # If gene is not in the master file, add it with TPM value for this sample
            echo -e "$gene_id\t$gene_name\t$tpm_value" >> "$MASTER_TSV"
        else
            # If gene is already in the master file, append the TPM value for this sample
            sed -i "/^$gene_id/s/$/\t$tpm_value/" "$MASTER_TSV"
        fi
    done < temp_file

    # Clean up temporary file
    rm temp_file
done

echo "StringTie processing, merging of GFF files and gene abundance files completed successfully."

# Define the input and output paths
VF_DB="/work_beegfs/sukem113/rUTI/Final/Virulence_Factor/VFDB_setB_nt.fas/VFDB_setB_nt.fas"
OUTPUT_VF="/work_beegfs/sukem113/rUTI/Final/Virulence_Factor/VFDB_setB_nt.fas/UTI89_only_virulence_factors.fas"

# Extract only lines with "Escherichia coli UTI89" entries
grep "^>.*Escherichia coli UTI89" "$VF_DB" > "$OUTPUT_VF"

echo "Filtered virulence factor database for Escherichia coli UTI89 created at $OUTPUT_VF"

# Define the input and output paths
FILTERED_VF_DB="/work_beegfs/sukem113/rUTI/Final/Virulence_Factor/VFDB_setB_nt.fas/UTI89_only_virulence_factors.fas"
OUTPUT_CSV="/work_beegfs/sukem113/rUTI/Final/Virulence_Factor/VFDB_setB_nt.fas/UTI89_gene_virulence_factors.csv"

# Drop the first and last columns
awk -F ' ' '{for(i=2; i<NF; i++) printf $i " "; print ""}' "$FILTERED_VF_DB" > "$OUTPUT_CSV"

echo "Output with first and last columns dropped saved to $OUTPUT_CSV"

The default line height has been increased for improved accessibility. You can choose to enable a more compact line height from the view settings menu.

‎Scripts/Main_Figures/Figure3.R
+45
Original file line number	Diff line number	Diff line change
@@ -0,0 +1,45 @@
library(tidyverse)
library(RColorBrewer)  # For accessing color palettes
library(scales)        # For formatting the axis labels
# Assuming gene_data is already loaded
# Set the path to your data and read in the data
file_path <- "consolidated_gene_data.csv"
gene_data <- read_csv(file_path)
# Filter out rows with FPKM == 0 to avoid -Inf in log scale
gene_data <- gene_data %>% filter(FPKM > 0)
# Define a custom, more saturated color palette
# Adjust these colors or add more based on the number of vf_category levels in your data
custom_palette <- c(
  "#D73027", "#FC8D59", "#FEE090", "#91BFDB", "#4575B4",
  "#1A9850", "#66BD63", "#A6D96A", "#D9EF8B", "#FFFFBF",
  "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F", "#9E0142"
)
# Create a scatter plot with adjusted color saturation
p <- ggplot(gene_data, aes(x = Sample, y = reorder(`Gene Name`, desc(`Gene Name`)), color = vf_category, size = FPKM)) +
  geom_point() +
  scale_size_continuous(range = c(1, 5)) +
  scale_color_manual(values = custom_palette) +
  theme_minimal() +
  labs(title = "",
       y = "Gene Name",
       x = "Sample",
       color = "Virulence Factor") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, face="bold", size=12),  # Bold and larger x-axis labels
        axis.text.y = element_text(face="bold", size=12),  # Bold and larger y-axis labels
        plot.title = element_text(hjust = 0.5, face="bold", size=14),  # Bold and larger plot title
        axis.title.x = element_text(face="bold", size=14),  # Bold and larger x-axis title
        axis.title.y = element_text(face="bold", size=12),  # Bold and larger y-axis title
        legend.title = element_text(face="bold", size=14),  # Bold and larger legend title
        legend.text = element_text(face="bold", size=12)) +  # Bold and larger legend text
  guides(color = guide_legend(override.aes = list(size = 5)))  # Adjust legend dot size for better visibility
# Display the plot
print(p)
# Save the plot
output_file_path <- "gene_expression_plot.png"
ggsave(filename = output_file_path, plot = p, width = 15, height = 10, units = "in")
