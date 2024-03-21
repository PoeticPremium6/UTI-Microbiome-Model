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
