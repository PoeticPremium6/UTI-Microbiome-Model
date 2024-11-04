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
output_file_path_genera <- "Figure2_genera.png"
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
  geom_bar(stat = "identity", color = "black", width = 1, alpha = 0) +  # Add black lines around bars for clarity
  labs(x = "Sample", y = "Relative Abundance (%) - Top 20 Species") +
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
output_file_path_species <- "Figure2_species.png"
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
ggsave("PCA_Samples_Composition_Improved.png", pca_samples_plot, width = 10, height = 8, dpi = 300)

#Is there a correlation between Enterobacter Diversity
#And PC1?
library(vegan)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(phyloseq)

# Filter physeq for Enterobacteriaceae family
entero_physeq <- subset_taxa(physeq, Family == "Enterobacteriaceae")

# Calculate Shannon diversity index for Enterobacteriaceae in each sample
entero_diversity <- estimate_richness(entero_physeq, measures = "Shannon")
colnames(entero_diversity) <- c("EnterobacteriaceaeDiversity")

# Merge Enterobacteriaceae diversity with PCA scores
# Ensure that your scores dataframe has "Sample" column matching rownames from entero_diversity
scores$Sample <- rownames(scores)
merged_df <- merge(scores, entero_diversity, by.x = "Sample", by.y = "row.names")

# Perform linear regression analysis
lm_result <- lm(EnterobacteriaceaeDiversity ~ PC1, data = merged_df)
summary_lm <- summary(lm_result)

# Extract R-squared value
r_squared <- summary_lm$r.squared

# Plot Enterobacteriaceae Diversity vs. PC1 Score with Sample Names and annotate with R-squared
p <- ggplot(merged_df, aes(x = PC1, y = EnterobacteriaceaeDiversity, label = Sample)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  geom_text_repel(size = 3, max.overlaps = Inf) +
  theme_minimal() +
  labs(title = "Enterobacteriaceae Diversity vs. PC1 Score",
       x = "PC1 Score",
       y = "Enterobacteriaceae Diversity") +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = Inf, y = -Inf, label = sprintf("R² = %.2f", r_squared), 
           hjust = 1, vjust = 0, size = 5, color = "blue")

# Display the plot
print(p)
# Assuming 'p' is the plot variable from the previous code
ggsave("Enterobacteriaceae_Diversity_vs_PC1.png", p, width = 10, height = 8, dpi = 300)

# Ensure you have patchwork installed
library(patchwork)

# Assuming pca_samples_plot and p are already defined 
#Define the layout to make the combined plot more square-like
combined_plot <- pca_samples_plot / p +
  plot_layout(heights = c(1, 1))

combined_plot <- combined_plot + plot_layout(guides = "collect") & 
  theme(aspect.ratio = 1)
# Display the combined plot
print(combined_plot)

# Save the combined plot to a file
ggsave("Combined_PCA_and_Diversity_Plot.png", combined_plot, width = 10, height = 16, dpi = 300)

                                          
# Load necessary library
library(ggplot2)

# Sample names
sample_names <- c("A01", "A02", "B01", "B02", "C01", "C02", "D01", "D02",
                  "E01", "E02", "F01", "F02", "G01", "H01", "H25361",
                  "H25362", "H25363", "H25364", "H25365")

# Shannon diversity values
shannon_div <- c(2.0695388, 0.7644332, 1.1369565, 1.6857164, 1.7876022,
                 0.9545996, 1.4798002, 1.3059761, 0.9599115, 1.1656087,
                 1.1563878, 0.7621552, 2.0830254, 1.7201654, 1.0412371,
                 0.9279537, 1.0937679, 0.9792370, 1.1679389)

# Combine into a data frame
shannon_div_df <- data.frame(
    Sample = sample_names,
    Shannon = shannon_div,
    stringsAsFactors = FALSE  
)

# Create the plot
plot <- ggplot(shannon_div_df, aes(x = Sample, y = Shannon)) +
    geom_bar(stat = "identity", fill = "purple") +
    labs(x = "Sample Names", y = "Alpha Diversity Measure (Shannon Index)") +
    theme(axis.title.x = element_text(face = "bold"),  # Bold x-axis title
          axis.title.y = element_text(face = "bold"),  # Bold y-axis title
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold"),  # Vertical x-axis ticks
          axis.text.y = element_text(face = "bold")) +  # Bold y-axis ticks
    theme_minimal()

# Save the plot
ggsave("/shannon_diversity_plot.png",
       plot, width = 10, height = 6)

