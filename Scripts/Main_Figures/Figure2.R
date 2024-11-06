#####Figure 2#####
# Define Paths to Data


# Load necessary libraries
library(phyloseq)
library(ggplot2)
library(dplyr)

# Custom color palette function to generate a diverse set of colors
generate_custom_palette <- function(n) {
  hues <- seq(15, 375, length.out = n + 1)
  hcl(h = hues, l = 40, c = 100)[1:n]
}

path_to_OTU <- "/rRNA_Counts_Clean.tsv"
path_to_TAX <- "/rRNA_Tax_Clean.tsv"

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

# Extract top genera
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
  labs(x = "Sample", y = "Relative Abundance (%)") +
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

# Save the plot for top  genera
output_file_path_genera <- "Figure2_genera.png"
ggsave(filename = output_file_path_genera, plot = p_genera, width = 10, height = 6)

# Display the plot
print(p_genera)

# Extract top  species
top_species <- rel_abundance_long %>%
  group_by(Species) %>%
  summarize(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance)) %>%
  top_n(20) %>%
  pull(Species)

# Filter data for top  species
filtered_rel_abundance_species <- rel_abundance_long %>%
  filter(Species %in% top_species)

# Generate a custom color palette using a diverse set of colors
color_palette_species <- generate_custom_palette(20)

# Plotting relative abundance for top  species with black lines between each
p_species <- ggplot(filtered_rel_abundance_species, aes(x = Sample, y = Abundance, fill = Species)) +
  geom_bar(stat = "identity") +
  geom_bar(stat = "identity", color = "black", width = 1, alpha = 0) +  # Add black lines around bars for clarity
  labs(x = "Sample", y = "Relative Abundance (%)") +
  scale_fill_manual(values = color_palette_species) +
  scale_y_continuous(breaks = seq(0, 100, by = 25), labels = paste0(seq(0, 100, by = 25), "%")) + # Define y-axis breaks and labels
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 12),  # Rotate x-axis labels
        axis.text.y = element_text(face = "bold", size = 12),  # Bold y-axis labels
        axis.title = element_text(face = "bold", size = 14),  # Bold and larger axis titles
        legend.position = "right",  # Legend on the right side
        legend.box = "vertical",  # Legend orientation
        legend.title = element_text(size = 10, face = "bold"),  # Legend title size and bold face
        legend.text = element_text(size = 12, face = "bold"),  # Legend text size and bold face
        panel.grid.major.x = element_blank(),  # Remove major grid lines for x-axis
        panel.grid.minor.x = element_blank(),  # Remove minor grid lines for x-axis
        panel.grid.major.y = element_line(color = "black", size = 0.2),  # Ensure major grid lines for y-axis are visible
        panel.grid.minor.y = element_blank()) +  # Remove minor grid lines for y-axis
  guides(fill = guide_legend(title = "Species", title.position = "top", title.theme = element_text(face = "bold")))  # Bold legend title

p_species  # Display the plot

# Save the plot for top 20 species
output_file_path_species <- "/Figure2_species.png"
ggsave(filename = output_file_path_species, plot = p_species, width = 10, height = 6)

# Display the plot
print(p_species)

#PCA of Species Composition
library(vegan)
library(ggplot2)

# Transform the OTU counts to relative abundances for PCA
physeq_relabund <- transform_sample_counts(physeq, function(x) x / sum(x))

# Perform PCA on the relative abundance data
otu_mat <- otu_table(physeq_relabund)
# Assuming 'otu_mat' contains your species composition data
pca_res <- prcomp(t(otu_mat), center = TRUE, scale. = TRUE)
# Calculate the percentage of variance explained by the PCs
variance_explained <- round(pca_res$sdev^2 / sum(pca_res$sdev^2) * 100, 2)

#Extract sample scores
scores <- as.data.frame(pca_res$x)
scores$Sample <- rownames(pca_res$x)

# Calculate the proportion of variance explained by each principal component
variance_explained <- pca_res$sdev^2 / sum(pca_res$sdev^2)
# Convert to percentages
variance_explained_percent <- round(variance_explained * 100, 2)

pca_samples_plot <- ggplot(scores, aes(x = PC1, y = PC2, label = Sample)) +
  geom_point() +
  geom_text_repel(aes(label = Sample), size = 5.5, 
                  box.padding = 1, point.padding = 0.3, 
                  max.overlaps = 50) +  # Adjusting max.overlaps for better label arrangement
  theme_minimal() +
  labs(title = "",
       x = paste("PC1 (", format(variance_explained_percent[1], nsmall = 2), "% variance explained)", sep=""),
       y = paste("PC2 (", format(variance_explained_percent[2], nsmall = 2), "% variance explained)", sep="")) +
  theme(legend.position = "none",
        plot.margin = unit(c(1, 1, 1, 1), "cm"),
        axis.title.x = element_text(size = 14, face = "bold"),  # Making x-axis title bold and size 12
        axis.title.y = element_text(size = 14, face = "bold")) +  # Making y-axis title bold and size 12
  xlim(min(scores$PC1) - 1, max(scores$PC1) + 1) +  # Adjusting limits for x-axis
  ylim(min(scores$PC2) - 1, max(scores$PC2) + 1)  # Adjusting limits for y-axis

# Display the plot
print(pca_samples_plot)

# Optionally, save the plot to a file
ggsave("/PCA_Samples_Composition_Improved.png", pca_samples_plot, width = 10, height = 8, dpi = 300)

library(ggplot2)
library(dplyr)

# Reverse the alphabetical order of species in the normalized abundance data
otu_path <- "/rRNA_Counts_Clean.tsv"
tax_path <- "/rRNA_Tax_Clean.tsv"

# Load and Preprocess Data
otu <- read.delim(otu_path, row.names = 1, check.names = FALSE, sep = "\t")
otu_matrix <- as.matrix(otu)
mode(otu_matrix) <- "numeric"  # Ensure numeric values in OTU table

# Taxonomy table
tax <- read.delim(tax_path, check.names = FALSE, sep = "\t", row.names = 1)
tax_matrix <- as.matrix(tax)

# Check and align sample names
sample_names <- colnames(otu_matrix)
sample_data <- data.frame(SampleID = sample_names, row.names = sample_names)

# Create phyloseq objects
otu_table_ps <- otu_table(otu_matrix, taxa_are_rows = TRUE)
tax_table_ps <- tax_table(tax_matrix)
sample_data_ps <- sample_data(sample_data)

# Create the Phyloseq object
physeq <- phyloseq(otu_table_ps, tax_table_ps, sample_data_ps)

# Normalize abundance data to 1000 per sample
otu_normalized <- sweep(otu_matrix, 2, colSums(otu_matrix) / 1000, FUN = "/")

# Replace OTUs with species names from the taxonomy table
species_names <- tax[rownames(otu_normalized), "Species"]

# Remove any rows with NA in species names
non_na_species <- !is.na(species_names)
otu_normalized <- otu_normalized[non_na_species, ]
species_names <- species_names[non_na_species]

# Update row names to species names
rownames(otu_normalized) <- species_names

# Convert to data frame for saving as .csv
otu_normalized_df <- as.data.frame(otu_normalized)

# Save the normalized abundance table
output_csv_path <- "/Species_Abundance_Per_Sample_Normalized.csv"
write.csv(otu_normalized_df, file = output_csv_path, row.names = TRUE)

# Display the first few rows to confirm
print(head(otu_normalized_df))

otu_normalized_df$Species <- factor(rownames(otu_normalized_df), 
                                    levels http://127.0.0.1:38295/graphics/df049fd8-f917-4f4a-ae55-bac6e87bdbd4.png= sort(unique(rownames(otu_normalized_df)), decreasing = TRUE))

# Reshape normalized data for plotting
# Convert to long format to use "Sample" and "Species" in ggplot
library(reshape2)
long_otu_normalized <- melt(otu_normalized_df, variable.name = "Sample", value.name = "Abundance")
long_otu_normalized$Species <- factor(long_otu_normalized$Species, levels = levels(otu_normalized_df$Species))

# Create the plot with purple tiles indicating presence and black dots sized by abundance
p_dotplot <- ggplot(long_otu_normalized, aes(x = Sample, y = Species)) +
  # Purple background tiles for presence of species
  geom_tile(data = subset(long_otu_normalized, Abundance > 0),
            fill = "purple", color = NA, alpha = 0.3, width = 0.9, height = 0.9) +
  # Black dots with size based on abundance
  geom_point(aes(size = Abundance), color = "black") +
  # Customize size scale for abundance with labels in legend
  scale_size_continuous(name = "Abundance", range = c(2, 8), 
                        breaks = c(200, 500, 1000), 
                        labels = c("200", "500", "1000")) +  
  # Custom labels for x, y, and size axes
  labs(x = "Sample", y = "Species", size = "Abundance") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 12),  # Rotate x-axis labels
    axis.text.y = element_text(face = "bold", size = 12),  # Bold y-axis labels
    axis.title = element_text(face = "bold", size = 14),  # Bold and larger axis titles
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    legend.position = "right",  # Position legend on right
    legend.title = element_text(face = "bold", size = 12),  # Bold legend title
    legend.text = element_text(size = 10)  # Legend text size
  ) +
  # Ensure only relevant legends are shown
  guides(fill = guide_legend(title = "Presence", override.aes = list(fill = "purple", alpha = 0.3)), 
         size = guide_legend(title = "Abundance"))

# Display the plot
print(p_dotplot)

# Save the plot to the specified path
output_file_path_dotplot <- "/Figure2_species_dotplot_New.png"
ggsave(filename = output_file_path_dotplot, plot = p_dotplot, width = 10, height = 8)

# Load necessary libraries
library(phyloseq)
library(ggplot2)

# Calculate Shannon Diversity Index
shannon_div_df <- estimate_richness(physeq, measures = "Shannon")
shannon_div_df$Sample <- rownames(shannon_div_df)

# Plot Shannon Diversity
plot <- ggplot(shannon_div_df, aes(x = Sample, y = Shannon)) +
  geom_bar(stat = "identity", fill = "purple") +
  labs(x = "Sample Names", y = "Alpha Diversity Measure (Shannon Index)") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(face = "bold"),  # Bold x-axis title
    axis.title.y = element_text(face = "bold"),  # Bold y-axis title
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold"),  # Rotate x-axis ticks
    axis.text.y = element_text(face = "bold")  # Bold y-axis ticks
  )

# Display the plot
print(plot)

# Save the plot
output_file_path_shannon <- "/shannon_diversity_plot.png"
ggsave(filename = output_file_path_shannon, plot = plot, width = 10, height = 6)
