#####Figure 3#####
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
pca_data$MarkerLabel <- ifelse(pca_data$IsContext, "CONTEXT", "CONTROL")
pca_data$PairKey <- gsub("(Control|Context)$", "", pca_data$Condition) # Key for pairing Control and Context

# Custom color palette
custom_palette <- c(
  "#D55E00", "#E69F00", "#F0E442", "#009E73", "#56B4E9",
  "#0072B2", "#CC79A7", "#999999", "#FF6666", "#77AC30",
  "#EDB120", "#7E2F8E", "#4DBEEE", "#A2142F", "#6666FF",
  "#FF00FF", "#00FFFF", "#CCFF00", "#FF9999"
)

# Prepare PCA Plot
p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_line(aes(group = PairKey), alpha = 0.5) +  # Draw lines between pairs
  geom_point(aes(shape = MarkerLabel), size = 3) +
  scale_shape_manual(values = c("CONTROL" = 16, "CONTEXT" = 17)) +
  scale_color_manual(values = custom_palette) +
  labs(x = "PC1", y = "PC2") +
  theme_minimal() +
  theme(axis.title.x = element_text(face = "bold"), axis.title.y = element_text(face = "bold"),
        legend.position = "right") +
  coord_fixed(ratio = 1.5)

# Bar Plot Data Preparation with Ordered Samples
bar_data <- data %>%
  pivot_longer(cols = -Metabolites, names_to = "Sample", values_to = "Value") %>%
  mutate(Group = gsub("(_.*)", "", Sample),
         Order = ifelse(grepl("Control", Sample), 1, 2)) %>%
  arrange(Group, Order) %>%
  group_by(Group, Sample) %>%
  summarise(Value = sum(Value, na.rm = TRUE), .groups = 'drop') %>%
  ungroup() %>%
  mutate(Sample = factor(Sample, levels = unique(Sample)))  # Ensuring the order is maintained

# Bar Plot with Vertical y-axis Labels and Ordered Samples
p_bar <- ggplot(bar_data, aes(x = Sample, y = Value, fill = Group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = custom_palette) +
  theme_minimal() +
  labs(x = "Sample", y = "Total Metabolic Flux / Biomass") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10, face = "bold"),
        axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),
        legend.position = "right")

# Combine PCA and Bar plots side by side without collecting legends
combined_plot <- p_pca + p_bar + plot_layout(ncol = 2)

# Save the Combined Plot
ggsave("Combined_PCA_Bar_Individual_Legends_Adjusted.png", combined_plot, width = 14, height = 7)
