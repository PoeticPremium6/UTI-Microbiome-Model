# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(patchwork)

# Load the data
data <- read.csv("Biomass.csv")

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
ggsave("Combined_PCA_Bar_Euclidean_Distance.png", combined_plot, width = 14, height = 10)
