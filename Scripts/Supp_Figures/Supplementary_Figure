###SUPPLEMENTARY FIGURES

#################################
#################################
###SUPP FIGURE
#Supp Figure 1
#Compare OTUs and TUs
library(ggplot2)
library(tibble)
library(dplyr)

# Sequencing technology information
seq_tech_data <- tibble::tribble(
  ~Name, ~Sequencing_Technology,
  "H25361", "Illumina MiSeq",
  "H25362", "Illumina MiSeq",
  "H25363", "Illumina MiSeq",
  "H25364", "Illumina MiSeq",
  "H25365", "Illumina MiSeq",
  "A01", "Illumina NovaSeq",
  "A02", "Illumina NovaSeq",
  "B01", "Illumina NovaSeq",
  "B02", "Illumina NovaSeq",
  "C01", "Illumina NovaSeq",
  "C02", "Illumina NovaSeq",
  "D01", "Illumina NovaSeq",
  "D02", "Illumina NovaSeq",
  "E01", "Illumina NovaSeq",
  "E02", "Illumina NovaSeq",
  "F01", "Illumina NovaSeq",
  "F02", "Illumina NovaSeq",
  "G01", "Illumina NovaSeq",
  "H01", "Illumina NovaSeq"
)

perform_pca_on_samples <- function(data_file, input_dir) {
  # Construct the full path to the data file
  full_path <- file.path(input_dir, data_file)
  
  # Read the data
  df <- read.csv(full_path, row.names = 1)
  
  # Transpose the dataframe so that samples are rows and features are columns
  df_t <- t(df)
  
  # Log-transform the data to improve normality
  df_log <- log1p(df_t)
  
  # Perform PCA
  pca_result <- prcomp(df_log, center = TRUE, scale. = TRUE)
  
  # Prepare a dataframe for ggplot
  pca_df <- as.data.frame(pca_result$x)
  pca_df$Sample <- row.names(pca_result$x)
  
  # Merge PCA data with sequencing technology information
  pca_df <- left_join(pca_df, seq_tech_data, by = c("Sample" = "Name"))
  
  # Plot PCA: PC1 vs PC2, color by Sequencing Technology
  ggplot_object <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Sequencing_Technology)) +
    geom_point(size = 4) +  # Increase the size of the points
    xlab(paste0("PC1 - ", round((pca_result$sdev[1]^2 / sum(pca_result$sdev^2)) * 100, 2), "% Variance")) +
    ylab(paste0("PC2 - ", round((pca_result$sdev[2]^2 / sum(pca_result$sdev^2)) * 100, 2), "% Variance")) +
    ggtitle(paste("PCA of Samples in", data_file)) +
    theme_minimal() +
    theme(
      legend.title = element_text(size = 12),  # Increase legend title size
      legend.text = element_text(size = 10),  # Increase legend text size
      plot.title = element_text(size = 14, face = "bold"),  # Bold and larger plot title
      plot.margin = unit(c(1, 1, 1, 1), "cm"),  # Adjust plot margins
      axis.text = element_text(size = 12),  # Larger axis text
      axis.title = element_text(size = 12)  # Larger axis title
    ) +
    scale_color_manual(values = c("Illumina MiSeq" = "blue", "Illumina NovaSeq" = "red")) +
    geom_text_repel(aes(label = Sample), size = 4, box.padding = 0.5)  # Use geom_text_repel for non-overlapping labels
  
  # Save the plot to the input directory
  ggsave(filename = file.path(input_dir, paste0(sub(".csv", "", data_file), "_PCA_plot.png")), plot = ggplot_object, width = 10, height = 8)
  
  print(ggplot_object)
}

input_dir <- "/UTI/"

# Perform PCA on OTU.csv and TU.csv
perform_pca_on_samples("OTU.csv", input_dir)
perform_pca_on_samples("TU.csv", input_dir)
#################################
#####Supp Figure 2#####
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

path_to_OTU <- "rRNA_Counts_Clean.tsv"
path_to_TAX <- "rRNA_Tax_Clean.tsv"

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

# Print Summary of the Phyloseq Object
print(physeq)

# Normalize the OTU table (relative abundance)
physeq_normalized <- transform_sample_counts(physeq, function(x) x / sum(x))

# Convert the normalized phyloseq object to a data frame
otu_normalized <- as.data.frame(otu_table(physeq_normalized))

# Apply a minimum threshold to define presence (e.g., 0.01), without setting values to 0
threshold <- 0.01

# Get strain names from the tax_table
strain_names <- as.character(tax_table(physeq)[, "Species"])

# Combine strain names with the normalized OTU table
otu_normalized_with_species <- cbind(Strain = strain_names, otu_normalized)

# Remove duplicate rows based on the 'Strain' column
otu_normalized_unique <- otu_normalized_with_species %>%
  distinct(Strain, .keep_all = TRUE)

# Adjust the filtering condition to be more inclusive
# Include species with at least one non-zero value above a smaller threshold (e.g., 0.0001)
inclusive_threshold <- 0.001
otu_normalized_filtered <- otu_normalized_unique %>%
  filter(rowSums(select(., -Strain) > inclusive_threshold) > 0)

# Save the normalized coverage to a CSV file
write.csv(otu_normalized_filtered, "normCoverage.csv", row.names = FALSE)

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
color_palette_genera <- c(
  "#4B0082", "#6A5ACD", "#483D8B", "#7B68EE", "#9370DB", 
  "#8A2BE2", "#9400D3", "#9932CC", "#BA55D3", "#DDA0DD", 
  "#EE82EE", "#8B008B", "#6B8E23", "#5F9EA0", "#4682B4",
  "#6495ED", "#00CED1", "#40E0D0", "#00BFFF", "#1E90FF"
)

# Create the plot
p_genera <- ggplot(filtered_rel_abundance, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  geom_bar(stat = "identity", color = "black", width = 1, alpha = 0) +  # Add black lines
  labs(x = "Sample", y = "Relative Abundance (%) - Top 20") +
  scale_fill_manual(values = color_palette_genera) +  # Apply the new color palette
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold"),  # Rotate x-axis labels
    axis.text.y = element_text(face = "bold"),  # Bold y-axis labels
    legend.position = "right",  # Legend on the right side
    legend.box = "vertical",  # Legend orientation
    legend.title = element_text(size = 10, face = "bold"),  # Legend title size and bold face
    legend.text = element_text(size = 8, face = "bold"),  # Legend text size and bold face
    panel.grid.major.x = element_blank(),  # Remove major grid lines
    panel.grid.minor.x = element_blank(),  # Remove minor grid lines
    panel.grid.major.y = element_line(color = "gray", size = 0.2),  # Add gray lines for major grid
    panel.grid.minor.y = element_blank()  # Remove minor grid lines
  ) +
  guides(fill = guide_legend(title = "Genus", title.position = "top", title.theme = element_text(face = "bold")))  # Bold legend title

# Save the plot for top 20 genera
output_file_path_genera <- "Supp2.png"
ggsave(filename = output_file_path_genera, plot = p_genera, width = 10, height = 6)

# Display the plot
print(p_genera)

#################################
#################################
#Supplementary Figure S3

# Load the required library
library(tidyverse)

# Read the non-context and context-specific model data
non_context_data <- read.csv("Combined_Metabolite_Uptake_Production_Summary.csv")
context_data <- read.csv("Combined_Metabolite_Uptake_Production_Summary.csv")

# Reshape non-context data to long format
non_context_long <- non_context_data %>%
  pivot_longer(cols = -Metabolite, names_to = "Sample", values_to = "NonContext_Amount")

# Reshape context data to long format
context_long <- context_data %>%
  pivot_longer(cols = -Metabolite, names_to = "Sample", values_to = "Context_Amount")

# Merge the long data by Metabolite and Sample
merged_data <- merge(non_context_long, context_long, by = c("Metabolite", "Sample"))

# Create a new column to distinguish the condition
merged_data$Condition <- ifelse(grepl("A01|A02|B01|B02", merged_data$Sample), "Control", "Context")

# Create a column for each sample/condition pair with clean column names
merged_data$Sample_Condition <- gsub("samp_", "", merged_data$Sample)  # Clean up the sample name
merged_data$Sample_Condition <- gsub("\\.mat", "", merged_data$Sample_Condition)  # Remove the .mat extension
merged_data$Sample_Condition <- gsub("microbiota_model_", "", merged_data$Sample_Condition)  # Clean prefix

# Separate NonContext and Context into distinct columns for better clarity
merged_data <- merged_data %>%
  mutate(NonContext_Column = paste(Sample_Condition, "NonContext", sep = "_"),
         Context_Column = paste(Sample_Condition, "Context", sep = "_"))

# Pivot wider, keeping NonContext and Context separate
data_noncontext <- merged_data %>%
  select(Metabolite, NonContext_Column, NonContext_Amount) %>%
  pivot_wider(names_from = NonContext_Column, values_from = NonContext_Amount)

data_context <- merged_data %>%
  select(Metabolite, Context_Column, Context_Amount) %>%
  pivot_wider(names_from = Context_Column, values_from = Context_Amount)

# Combine the NonContext and Context data
final_data <- full_join(data_noncontext, data_context, by = "Metabolite")

# Print the reshaped data
print(final_data)

# Ensure the target directory exists
output_directory <- "Supp/"

# Specify the output file name (you can adjust the name as needed)
output_file <- paste0(output_directory, "reshaped_data.csv")

# Save the reshaped data to the specified directory
#Adjusted file head to 'Control'
write.csv(final_data, file = output_file, row.names = FALSE)

# Confirm that the data was saved
cat("Data saved to:", output_file, "\n")
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(patchwork)

# Load the data
data <- read.csv("/Supp/reshaped_data.csv")

# Transpose and prepare data for PCA
data_for_pca <- t(data[, -1])
nonconst_cols <- apply(data_for_pca, 2, var) > 0
data_for_pca <- data_for_pca[, nonconst_cols]

# Perform PCA
pca_result <- prcomp(data_for_pca, scale. = TRUE)
pca_data <- as.data.frame(pca_result$x[, 1:2])
pca_data$Condition <- rownames(pca_data)
pca_data$Group <- gsub("(\\D+)(\\d+).*", "\\1\\2", pca_data$Condition)
pca_data$IsContext <- grepl("Context", pca_data$Condition)

# Adjust MarkerLabel to "Context" and "non-Context"
pca_data$MarkerLabel <- ifelse(pca_data$IsContext, "Context", "non-Context")

pca_data$PairKey <- gsub("(Control|Context)$", "", pca_data$Condition) # Key for pairing Control and Context

# Custom color palette
custom_palette <- c(
  "#D55E00", "#E69F00", "#F0E442", "#009E73", "#56B4E9",
  "#0072B2", "#CC79A7", "#999999", "#FF6666", "#77AC30",
  "#EDB120", "#7E2F8E", "#4DBEEE", "#A2142F", "#6666FF",
  "#FF00FF", "#00FFFF", "#CCFF00", "#FF9999"
)

# Calculate the proportion of variance explained by each principal component
variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
variance_explained_percent <- variance_explained * 100

# Update the PCA plot to include variance explained in the axis titles
p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_line(aes(group = PairKey), alpha = 0.5) +  # Draw lines between pairs
  geom_point(aes(shape = MarkerLabel), size = 4) +
  scale_shape_manual(values = c("non-Context" = 16, "Context" = 17)) +
  scale_color_manual(values = custom_palette) +
  labs(x = paste("PC1 (", sprintf("%.2f%%", variance_explained_percent[1]), " variance explained)", sep=""),
       y = paste("PC2 (", sprintf("%.2f%%", variance_explained_percent[2]), " variance explained)", sep="")) +
  theme_minimal() +
  theme(axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        legend.title = element_text(size = 12, face = "bold"),
        legend.position = "right",
        legend.text = element_text(face = "bold", size = 12)) +
  coord_fixed(ratio = 1.5)

# Display the updated PCA plot
print(p_pca)

# Bar Plot Data Preparation with Ordered Samples
library(ggplot2)
library(dplyr)
library(tidyr)

# Assuming 'data' has been loaded
# Prepare the bar_data with proper group naming
bar_data <- data %>%
  pivot_longer(cols = -Metabolites, names_to = "Sample", values_to = "Value") %>%
  mutate(
    Group = ifelse(grepl("Control", Sample), "non-Context", "Context"), # Rename groups to "context" and "non-context"
    SampleID = gsub("(Control|Context)", "", Sample),
    SampleID = gsub("_$", "", SampleID) # Remove trailing underscore
  ) %>%
  group_by(Group, SampleID) %>%
  summarise(TotalFlux = sum(Value), .groups = 'drop') %>%
  ungroup()

# Ensure the SampleID factor levels are set correctly
bar_data$SampleID <- factor(bar_data$SampleID, levels = unique(bar_data$SampleID))

# Use your existing custom_palette for coloring
custom_palette <- c(
  "#D55E00", "#E69F00", "#F0E442", "#009E73", "#56B4E9",
  "#0072B2", "#CC79A7", "#999999", "#FF6666", "#77AC30",
  "#EDB120", "#7E2F8E", "#4DBEEE", "#A2142F", "#6666FF",
  "#FF00FF", "#00FFFF", "#CCFF00", "#FF9999"
)

# Adjusted Bar Plot for Stacked Visualization with updated group names
p_bar_stacked_custom_palette <- ggplot(bar_data, aes(x = Group, y = TotalFlux, fill = SampleID)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = custom_palette[1:length(unique(bar_data$SampleID))]) + # Apply custom color palette with dynamic length based on SampleID
  theme_minimal() +
  labs(x = "", y = "Total Metabolic Flux / Biomass", fill = "Sample ID") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 12, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"), 
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.position = "right")

# Display the plot
print(p_bar_stacked_custom_palette)



# Combine PCA and Bar plots side by side without collecting legends
combined_plot <- p_pca + p_bar_stacked_custom_palette + plot_layout(ncol = 2)

# Save the Combined Plot
ggsave("Combined_PCA_Bar.png", combined_plot, width = 14, height = 7)

# Assuming pca_data is already prepared
library(dplyr)

library(tidyr)
library(dplyr)

# Assuming 'data' is your original dataframe
data_long <- pivot_longer(data, cols = -Metabolites, 
                          names_to = "Sample_Condition", 
                          values_to = "Flux")

# Separate Sample and Condition
data_long <- data_long %>%
  separate(Sample_Condition, into = c("Sample", "Condition"), sep = "_")

# Calculate total flux for each condition within each sample
total_flux <- data_long %>%
  group_by(Sample, Condition) %>%
  summarise(TotalFlux = sum(Flux), .groups = 'drop')

# Ensure total_flux is in the correct format
total_flux$Sample <- factor(total_flux$Sample)
total_flux$Condition <- factor(total_flux$Condition)

# Pivot wider to have Control and Context side-by-side
flux_wide <- pivot_wider(total_flux, names_from = Condition, values_from = TotalFlux)

# Calculate Euclidean distance (simple subtraction in this context as we're dealing with a single dimension)
flux_wide$EuclideanDistance <- abs(flux_wide$Control - flux_wide$Context)

library(ggplot2)

# Ensure there are enough colors for all samples
sample_colors <- setNames(custom_palette[1:length(unique(flux_wide$Sample))], levels(flux_wide$Sample))

ggplot(flux_wide, aes(x = Sample, y = EuclideanDistance, fill = Sample)) +
  geom_col() +
  scale_fill_manual(values = sample_colors) + # Apply the custom colors
  theme_minimal() +
  labs(title = "Euclidean Distance Between Control and Context Conditions",
       x = "Sample", y = "Euclidean Distance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "none") # Hide the legend if not needed

library(patchwork)

# Adjust the layout of individual plots if necessary
# For box-shaped plots at the top, you might consider adjusting their aspect ratio or size to fit your description
p_pca_design <- p_pca + plot_layout(heights = c(1))
p_bar_stacked_design <- p_bar_stacked_custom_palette + plot_layout(heights = c(1))

# For the Euclidean Distance plot, ensure it's created as per your description
p_euclidean_distance <- ggplot(flux_wide, aes(x = Sample, y = EuclideanDistance, fill = Sample)) +
  geom_col() +
  scale_fill_manual(values = sample_colors) +
  theme_minimal() +
  labs(title = "", x = "Sample", y = "Euclidean Distance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14, face = "bold"),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        legend.position = "none")


# Combine plots with specified layout
# Use `plot_layout()` to define the number of columns and the relative heights of the plots
combined_plot <- (p_pca_design | p_bar_stacked_design) / 
  p_euclidean_distance + 
  plot_layout(ncol = 1, heights = c(1, 0.5))

# Display the combined plot
print(combined_plot)

# Save the Combined Plot with appropriate dimensions
ggsave("Supp3.png", combined_plot, width = 14, height = 12)

