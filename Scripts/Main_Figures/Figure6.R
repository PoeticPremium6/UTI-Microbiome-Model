# Load necessary libraries
library(ggplot2)
library(reshape2)  # For melt function
library(viridis)   # For color scale
library(dplyr)     # For data manipulation
library(readr)     # For read_csv

# Read data
data <- read.csv("PA.csv")

# Filter to keep only metabolites showing both presence (1) and absence (0) across samples
data_filtered <- data %>%
  mutate(across(-Metabolite, as.numeric)) %>%
  filter(apply(., 1, function(x) all(c(0, 1) %in% x[-1])))

# Transform the filtered data into a long format
data_long <- melt(data_filtered, id.vars = "Metabolite")

# Reverse the levels on the y-axis
data_long$variable <- factor(data_long$variable, levels = rev(levels(data_long$variable)))

# Plotting with visual enhancements
p <- ggplot(data_long, aes(x = Metabolite, y = variable)) + 
  geom_tile(aes(fill = factor(value)), color = "white") + 
  scale_fill_viridis_d(name = "Presence", begin = 0.2, end = 0.8,
                       labels = c("0" = "Absent", "1" = "Present")) + 
  labs(x = "Metabolite", y = "Sample") + 
  theme_minimal(base_size = 14) + 
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 16, face = "bold"),  # Increase size and make bold for x-axis labels
    axis.text.y = element_text(size = 16, face = "bold"),  # Increase size and make bold for y-axis labels
    axis.title.x = element_text(size = 18, face = "bold"),  # Increase size and make bold for x-axis title
    axis.title.y = element_text(size = 18, face = "bold"),  # Increase size and make bold for y-axis title
    legend.position = "right",
    legend.title = element_text(size = 16, face = "bold"),  # Increase size and make bold for legend title
    legend.text = element_text(size = 14, face = "bold"),  # Increase size and make bold for legend text
    plot.margin = margin(1, 1, 1, 1, "cm")
  ) +
  coord_fixed(ratio = 1) +
  geom_hline(yintercept = seq(2.5, length(unique(data_long$variable)), by = 2), color = "black", size = 1.5)  # Adjusted line placement

# Display the plot
print(p)



# Save the plot in the specified directory with improved quality
ggsave("metabolite_heatmap_differences_filtered.png", p, width = 16, height = 12, units = "in")

# Load necessary libraries
library(ggplot2)
library(reshape2)  # For melt function
library(viridis)   # For color scale
library(dplyr)     # For data manipulation
library(readr)     # For read_csv
library(grid)      # For arrow type in geom_segment
library(ggrepel)   # For non-overlapping text labels

# Filter to keep only metabolites showing both presence (1) and absence (0) across samples
data_filtered <- data %>%
  mutate(across(-Metabolite, as.numeric)) %>%
  filter(apply(., 1, function(x) all(c(0, 1) %in% x[-1])))

# Preparing data for PCA: Exclude 'Metabolite' column and transpose
data_for_pca <- data_filtered[,-1]  # Excluding the Metabolite column for PCA
rownames(data_for_pca) <- data_filtered$Metabolite

# Transpose data for PCA
data_t <- t(data_for_pca)

# Convert data to numeric and normalize
data_norm <- scale(apply(data_t, 2, as.numeric))

# Perform PCA
pca_res <- prcomp(data_norm, center = TRUE, scale. = TRUE)

# Extract the proportion of variance explained by the first two principal components
variance_explained <- pca_res$sdev^2 / sum(pca_res$sdev^2)
pc1_explained <- variance_explained[1]
pc2_explained <- variance_explained[2]

# Convert loadings to a dataframe and ensure metabolite names are included
loadings_df <- as.data.frame(pca_res$rotation)
loadings_df$Metabolite <- rownames(pca_res$rotation)

# Plotting the loadings with metabolite names using ggrepel to avoid overlap
# Improved loading plot
loading_plot_improved <- ggplot(loadings_df, aes(x = PC1, y = PC2, label = Metabolite)) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(type = "closed", length = unit(0.15, "inches")), color = "red", size = 0.5) +
  geom_text_repel(aes(label = Metabolite), 
                  size = 4.5,  # Increase the text size
                  box.padding = unit(0.2, "lines"), 
                  point.padding = unit(0.3, "lines"), 
                  min.segment.length = unit(0.1, "inches"),
                  max.overlaps = 100) +  # Allow more overlaps to ensure all labels are printed
  theme_minimal() +
  labs(title = "Loadings (Metabolite Contributions)", 
       x = paste("Principal Component 1 (", sprintf("%.2f%%", pc1_explained * 100), " variance explained)", sep=""),
       y = paste("Principal Component 2 (", sprintf("%.2f%%", pc2_explained * 100), " variance explained)", sep="")) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.position = "none")  # Remove legend if unnecessary

# View the improved loading plot
print(loading_plot_improved)

# Save the biplot
ggsave("PCA_biplot.png", loading_plot_improved, width = 10, height = 8, dpi = 300)


# Get a summary of presence (1) for all samples
all_presence_summary <- data_filtered %>%
  pivot_longer(-Metabolite, names_to = "Sample", values_to = "Presence") %>%
  group_by(Metabolite) %>%
  summarise(Total_Presence = sum(Presence == 1, na.rm = TRUE)) %>%
  ungroup()

# Get presence summary specifically for G1 and E2
selected_samples_summary <- data_filtered %>%
  select(Metabolite, `G1_Control`, `G1_Context`, `E2_Control`, `E2_Context`) %>%
  pivot_longer(-Metabolite, names_to = "Sample", values_to = "Presence") %>%
  group_by(Metabolite) %>%
  summarise(Selected_Presence = sum(Presence == 1, na.rm = TRUE)) %>%
  ungroup()

# Find metabolites with unique presence in G1 and E2 by comparing to all samples
unique_metabolites <- selected_samples_summary %>%
  filter(Selected_Presence > 0) %>%
  left_join(all_presence_summary, by = "Metabolite") %>%
  filter(Total_Presence == Selected_Presence) %>%
  select(Metabolite)

# Print unique metabolites for G1 and E2
print(unique_metabolites)
