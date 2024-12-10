# Load necessary libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(stats)
library(readxl)

data_path <- "Subsystems_Merged.xlsx"
data <- read_excel(data_path)

# Reshape the data to long format for analysis
data_long <- data %>%
  pivot_longer(cols = -Subsystems, names_to = c("Sample_Condition"), values_to = "Activity") %>%
  separate(Sample_Condition, into = c("Sample", "Condition"), sep = "_")

# View the reshaped data
head(data_long)

# Perform a paired t-test (or Wilcoxon test if the data is not normally distributed)
test_results <- data_long %>%
  group_by(Subsystems) %>%
  summarise(
    p_value = t.test(Activity ~ Condition, paired = TRUE)$p.value
  )

# Adjust p-values for multiple comparisons using Benjamini-Hochberg
test_results$adj_p_value <- p.adjust(test_results$p_value, method = "BH")

# View the results with adjusted p-values
head(test_results)
# Filter for significant results (adjusted p-value < 0.05)
# Filter for significant results (adjusted p-value < 0.05)
significant_results <- test_results %>%
  filter(adj_p_value < 0.05)

significant_results <- significant_results %>%
  mutate(log_adj_p_value = -log10(adj_p_value))

# Create a plot of the significant subsystems only with log-transformed adjusted p-values
bar_plot <- ggplot(significant_results, aes(x = reorder(Subsystems, adj_p_value), y = adj_p_value)) +
  geom_bar(stat = "identity", fill = "purple") +  # Bar plot with purple color
  coord_flip() +  # Flip the coordinates for better readability
  labs(
    title = "", 
    x = "Subsystem", 
    y = "-Log10(Adjusted P-value)"
  ) +  # Use log-transformed label
  theme_minimal() +  # Minimal theme for the plot
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  # Bold x-axis ticks
    axis.text.y = element_text(face = "bold", size = 10),  # Increased font size and bold y-axis ticks
    axis.title.x = element_text(face = "bold", size = 12),  # Bold x-axis title
    axis.title.y = element_text(face = "bold", size = 12),  # Bold y-axis title
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)  # Centered bold title
  )

# Save the plot with tall dimensions
ggsave(
  filename = "Tall_Bar_Plot.png",
  plot = bar_plot,
  width = 6,  # Narrower width
  height = 10,  # Taller height
  dpi = 300
)

# Display the plot
print(bar_plot)
#########
# Load necessary libraries
library(dplyr)
library(tidyr)
library(pheatmap)
library(readxl)
library(grid)

# Step 0: Load the dataset
data_path <- "Subsystems_Merged.xlsx"
data <- read_excel(data_path)

# Step 1: Reshape the data to long format
data_long <- data %>%
  pivot_longer(cols = -Subsystems, names_to = c("Sample_Condition"), values_to = "Activity") %>%
  separate(Sample_Condition, into = c("Sample", "Condition"), sep = "_")

# Step 2: Reshape the data into wide format
data_wide <- data_long %>%
  pivot_wider(names_from = Condition, values_from = Activity) %>%
  mutate(fold_change = Context / NonContext) %>%
  select(Subsystems, Sample, NonContext, Context, fold_change)

# Step 3: Handle missing or infinite values
data_wide_clean <- data_wide %>%
  filter(!is.na(fold_change) & !is.infinite(fold_change))

# Step 4: Filter subsystems with low variance across samples
low_variance_threshold <- 0.30  # Threshold for variance filtering
high_value_subsystems <- data_wide_clean %>%
  group_by(Subsystems) %>%
  summarize(variance = var(fold_change, na.rm = TRUE)) %>%
  filter(variance > low_variance_threshold) %>%
  pull(Subsystems)

# Step 5: Keep only high-value subsystems
data_wide_filtered <- data_wide_clean %>%
  filter(Subsystems %in% high_value_subsystems)

# Step 6: Create a matrix for the heatmap
data_matrix <- data_wide_filtered %>%
  select(Sample, Subsystems, fold_change) %>%
  spread(key = Subsystems, value = fold_change)

# Ensure sample names are correctly set as rownames
rownames(data_matrix) <- data_matrix$Sample
data_matrix <- data_matrix %>% select(-Sample)  # Remove Sample column after setting rownames

# Step 7: Define an improved color palette
improved_palette <- colorRampPalette(c("green", "white", "purple"))(100)

library(grid)
# Step 8: Generate the heatmap with swapped x and y axes, keeping sample names on the y-axis
data_matrix_transposed <- t(data_matrix)  # Transpose the matrix
colnames(data_matrix_transposed) <- rownames(data_matrix)  # Keep sample names as column names

# Step 9: Generate the heatmap with bold axis labels and ticks
heatmap_plot <- pheatmap(
  data_matrix_transposed,  # Transposed matrix
  scale = "row",  # Normalize rows to focus on relative differences
  clustering_distance_rows = "euclidean",  # Hierarchical clustering for rows (subsystems)
  clustering_distance_cols = "euclidean",  # Hierarchical clustering for columns (samples)
  clustering_method = "complete",  # Clustering method
  color = improved_palette,  # Enhanced color palette
  fontsize = 10,  # General font size
  fontsize_row = 22,  # Font size for subsystem labels (now y-axis)
  fontsize_col = 12,  # Font size for sample labels (now x-axis)
  legend_breaks = seq(-2, 2, by = 0.5),  # Adjust breaks for color legend
  legend_labels = seq(-2, 2, by = 0.5),  # Corresponding labels for the color scale
  legend_title = "Fold Change (Log Scale)",  # Legend title
  legend_side = "left",  # Move legend to the left side
  angle_col = 45,  # Rotate subsystem names for better readability
  # Bold and larger font size for axis labels
  axis.text.x = grid::gpar(fontsize = 14, fontface = "bold"),  # Bold and larger font for x-axis ticks
  axis.text.y = grid::gpar(fontsize = 14, fontface = "bold"),  # Bold and larger font for y-axis ticks
  # Adjust color legend to match the color range and align it properly
  color_fun = colorRampPalette(c("green", "white", "purple")),  # Ensure colors match the range
  # Increase resolution and better positioning of color bar labels
  display_numbers = F  # Optionally, remove numbers inside heatmap cells for clarity
)

# Print the heatmap plot to verify appearance
heatmap_plot

# Print the heatmap plot to verify appearance
heatmap_plot

# Specify the file path where you want to save the figure
output_path <- "heatmap.png"

# Save the heatmap as a PNG file
png(output_path, width = 1000, height = 1000)  # Set dimensions as needed
heatmap_plot  # Plot the heatmap to the PNG file
dev.off()  # Close the device to save the file


###############
###############
# Load necessary libraries
library(tidyverse)
library(readxl)
library(ggpubr)

# Load the data
data_path <- "Subsystems_Merged.xlsx"
data <- read_excel(data_path)

# Reshape the data to long format
data_long <- data %>%
  pivot_longer(cols = -Subsystems, names_to = c("Sample_Condition"), values_to = "Activity") %>%
  separate(Sample_Condition, into = c("Sample", "Condition"), sep = "_") %>%
  mutate(
    Condition = str_trim(tolower(Condition))  # Standardize condition names
  )

# 1. Create a summary of activity data for each Lactobacillus group
summary_stats <- data_long %>%
  filter(Community %in% c("Lactobacillus Diverse", "Lactobacillus Absent", "Lactobacillus Single")) %>%
  group_by(Community, Condition) %>%
  summarise(
    mean_activity = mean(Activity, na.rm = TRUE),
    sd_activity = sd(Activity, na.rm = TRUE),
    min_activity = min(Activity, na.rm = TRUE),
    max_activity = max(Activity, na.rm = TRUE),
    n = n(),  # Sample size for each group
    .groups = 'drop'
  )

# Print summary statistics for each group and condition
print("Summary Statistics for Lactobacillus Groups:")
print(summary_stats)

# 2. Perform ANOVA to test for differences in activity between Lactobacillus groups
anova_results <- data_long %>%
  filter(Community %in% c("Lactobacillus Diverse", "Lactobacillus Absent", "Lactobacillus Single")) %>%
  aov(Activity ~ Community + Condition + Community:Condition, data = .)

# Summary of ANOVA results
anova_summary <- summary(anova_results)
print("ANOVA Results for Activity Across Lactobacillus Groups:")
print(anova_summary)

# 3. Perform pairwise comparisons using Tukey's HSD test (post-hoc test after ANOVA)
tukey_results <- TukeyHSD(anova_results)

# Print Tukey's HSD post-hoc test results
print("Tukey's HSD Post-hoc Test Results for Lactobacillus Groups:")
print(tukey_results)

# Define the communities
#Diverse communities are more than 3 species
Lactobacillus_diverse <- c("A01", "B02", "D01", "D02", "G01")
Lactobacillus_absent <- c("C02", "E01", "E02", "F01", "F02", "H01", "H25362", "H25363", "H25364", "H25365")
Lactobacillus_single <- c("A02", "B01", "C01", "H25361")

# Assign community type to each sample
data_long <- data_long %>%
  mutate(
    Community = case_when(
      Sample %in% Lactobacillus_diverse ~ "Lactobacillus Diverse",
      Sample %in% Lactobacillus_absent ~ "Lactobacillus Absent",
      Sample %in% Lactobacillus_single ~ "Lactobacillus Single",
      TRUE ~ "Other"
    )
  )

# Split the data into 'Context' and 'NonContext' datasets
data_context <- data_long %>% filter(Condition == "context")
data_noncontext <- data_long %>% filter(Condition == "noncontext")

# Custom function to create the dot plots with enhanced statistical test and comparisons
create_dot_plot <- function(data, condition_label) {
  ggplot(data, aes(x = Community, y = Activity, fill = Community)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.15, size = 2.5, aes(color = Community)) +
    scale_fill_manual(values = c("#4B0082", "#800080", "#D8BFD8")) +  # Deep purple to light purple palette
    scale_color_manual(values = c("#4B0082", "#800080", "#D8BFD8")) +  # Matching color for jitter points
    labs(
      title = paste("Metabolic Activity (", condition_label, ")", sep = ""),
      x = "Community Type",
      y = "Activity"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  # Centered bold title
      axis.text.x = element_text(size = 10, face = "bold"),  # Bold and larger x-axis ticks
      axis.text.y = element_text(size = 12, face = "bold"),  # Bold and larger y-axis ticks
      axis.title.y = element_text(size = 12, face = "bold"),  # Bold y-axis title
      axis.title.x = element_text(size = 12, face = "bold"),  # Bold x-axis title
      legend.position = "none"  # Remove legend for clarity
    ) +
    # Perform and display statistical tests (t-tests)
    stat_compare_means(
      comparisons = list(
        c("Lactobacillus Diverse", "Lactobacillus Absent"),
        c("Lactobacillus Diverse", "Lactobacillus Single"),
        c("Lactobacillus Absent", "Lactobacillus Single")
      ),
      method = "t.test",
      label = "p.signif",  # Show significance level (e.g., *, **, ***) 
      size = 4,
      position = position_dodge(width = 0.5),  # Adjust position of comparison labels to prevent overlap
      label.x.npc = "center",  # Center comparison labels horizontally
      label.y = 1.2,  # Adjust vertical position to avoid overlap
      angle = 45,  # Rotate labels by 45 degrees
      tip_length = 0.03  # Adjust the tip length of the lines connecting groups
    )
}

# Create the plots for context and noncontext
context_plot <- create_dot_plot(data_context, "Context")
noncontext_plot <- create_dot_plot(data_noncontext, "NonContext")

# Combine the plots into a single figure (side by side)
combined_plot <- ggarrange(context_plot, noncontext_plot, ncol = 2, labels = c("A", "B"))

# Save the combined plot with the specified size
ggsave("Combined_Plot_Context_NonContext.png", 
       combined_plot, width = 12, height = 6, dpi = 300)

# Display the combined plot
print(combined_plot)
