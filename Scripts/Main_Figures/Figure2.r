#####Figure 2#####
# Load necessary libraries
library(phyloseq)
library(ggplot2)
library(dplyr)

# Custom color palette function to generate a diverse set of colors
generate_custom_palette <- function(n) {
  hues <- seq(15, 375, length.out = n + 1)
  hcl(h = hues, l = 40, c = 100)[1:n]
}

path_to_OTU <- "/Users/josspa/UTI/Taxonomy/Merged_Counts.tsv"
path_to_TAX <- "/Users/josspa/UTI/Taxonomy/Merged_Tax.tsv"

# Import the OTU data
OTU <- read.delim(path_to_OTU, header = TRUE, sep = "\t", check.names = FALSE)

# Handle duplicate or missing OTU IDs
OTU$OTU[OTU$OTU == ""] <- paste0("EmptyID_", seq_len(sum(OTU$OTU == "")))
dup_ids <- which(duplicated(OTU$OTU))
OTU$OTU[dup_ids] <- paste0(OTU$OTU[dup_ids], "_dup", seq_len(length(dup_ids)))

# Set OTU IDs as row names and convert to matrix
rownames(OTU) <- OTU$OTU
OTU$OTU <- NULL  
OTU <- otu_table(as.matrix(OTU), taxa_are_rows = TRUE)

# Import the Taxonomy data
TAX <- read.delim(path_to_TAX, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
TAX <- tax_table(as.matrix(TAX))

# Find common taxa between OTU and TAX tables
common_taxa <- intersect(rownames(OTU), rownames(TAX))

# Subset OTU and TAX tables to only contain common taxa
OTU <- OTU[common_taxa, ]
TAX <- TAX[common_taxa, ]

# Ensure that OTU and Taxonomy tables are compatible (they have matching row names)
if(!all(rownames(OTU) == rownames(TAX))) stop("Mismatch between OTU and TAX tables")

# Manually creating sample data
sample_names <- c("A01", "A02", "B01", "B02", "C01", "C02", "D01", "D02",
                  "E01", "E02", "F01", "F02", "G01", "H01", "H25361",
                  "H25362", "H25363", "H25364", "H25365")

SAMPLES <- data.frame(SampleID = sample_names, Group = rep("A", length(sample_names))) # Group column is arbitrary
row.names(SAMPLES) <- sample_names
SAMPLES <- sample_data(SAMPLES)

# Create phyloseq object
physeq <- phyloseq(OTU, TAX, SAMPLES)

# Calculate relative abundance
rel_abundance <- transform_sample_counts(physeq, function(x) x / sum(x))

# Convert to long format for plotting
rel_abundance_long <- psmelt(rel_abundance)

# Extract top 20 genera
top_genera <- rel_abundance_long %>%
  group_by(Genus) %>%
  summarize(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance)) %>%
  top_n(20) %>%
  pull(Genus)

# Filter data for top 20 genera
filtered_rel_abundance <- rel_abundance_long %>%
  filter(Genus %in% top_genera)

# Generate a custom color palette using a diverse set of colors
color_palette_genera <- generate_custom_palette(20)

# Plotting relative abundance for top 20 genera with black lines between each
p_genera <- ggplot(filtered_rel_abundance, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  geom_bar(stat = "identity", color = "black", width = 1, alpha = 0) +  # Add black lines
  labs(x = "Sample", y = "Relative Abundance (%) - Top 20") +
  scale_fill_manual(values = color_palette_genera) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold"),  # Rotate x-axis labels
        axis.text.y = element_text(face = "bold"),  # Bold y-axis labels
        legend.position = "right",  # Legend on the right side
        legend.box = "vertical",  # Legend orientation
        legend.title = element_text(size = 10, face = "bold"),  # Legend title size and bold face
        legend.text = element_text(size = 8, face = "bold"),  # Legend text size and bold face
        panel.grid.major.x = element_blank(),  # Remove major grid lines
        panel.grid.minor.x = element_blank(),  # Remove minor grid lines
        panel.grid.major.y = element_line(color = "gray", size = 0.2),  # Add gray lines for major grid
        panel.grid.minor.y = element_blank()) +  # Remove minor grid lines
  guides(fill = guide_legend(title = "Genus", title.position = "top", title.theme = element_text(face = "bold")))  # Bold legend title

# Save the plot for top 20 genera
output_file_path_genera <- "/Users/josspa/UTI/Figures/Figure2/Figure2_genera.png"
ggsave(filename = output_file_path_genera, plot = p_genera, width = 10, height = 6)

# Display the plot
print(p_genera)

# Extract top 20 species
top_species <- rel_abundance_long %>%
  group_by(Species) %>%
  summarize(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance)) %>%
  top_n(20) %>%
  pull(Species)

# Filter data for top 20 species
filtered_rel_abundance_species <- rel_abundance_long %>%
  filter(Species %in% top_species)

# Generate a custom color palette using a diverse set of colors
color_palette_species <- generate_custom_palette(20)

# Plotting relative abundance for top 20 species with black lines between each
p_species <- ggplot(filtered_rel_abundance_species, aes(x = Sample, y = Abundance, fill = Species)) +
  geom_bar(stat = "identity") +
  geom_bar(stat = "identity", color = "black", width = 1, alpha = 0) +  # Add black lines
  labs(x = "Sample", y = "Relative Abundance (%) - Top 20") +
  scale_fill_manual(values = color_palette_species) +
  scale_y_continuous(breaks = seq(0, 100, by = 25), labels = paste0(seq(0, 100, by = 25), "%")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold"),  # Rotate x-axis labels
        axis.text.y = element_text(face = "bold"),  # Bold y-axis labels
        legend.position = "right",  # Legend on the right side
        legend.box = "vertical",  # Legend orientation
        legend.title = element_text(size = 10, face = "bold"),  # Legend title size and bold face
        legend.text = element_text(size = 8, face = "bold"),  # Legend text size and bold face
        panel.grid.major.x = element_blank(),  # Remove major grid lines
        panel.grid.minor.x = element_blank(),  # Remove minor grid lines
        panel.grid.major.y = element_line(color = "gray", size = 0.2),  # Add gray lines for major grid
        panel.grid.minor.y = element_blank()) +  # Remove minor grid lines
  guides(fill = guide_legend(title = "Species", title.position = "top", title.theme = element_text(face = "bold")))  # Bold legend title

# Save the plot for top 20 species
output_file_path_species <- "/Users/josspa/UTI/Figures/Figure2/Figure2_species.png"
ggsave(filename = output_file_path_species, plot = p_species, width = 10, height = 6)

# Display the plot
print(p_species)

