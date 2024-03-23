library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(scales) 

# Read the data
data <- read.csv("Class_Clustering.csv")

# Keep only SuperClass and sample columns
data_sum <- data %>%
  select(-Class, -SubClass, -Metabolite) %>%
  pivot_longer(cols = starts_with("A") | starts_with("B") | starts_with("C") | starts_with("D") | starts_with("E") | starts_with("F") | starts_with("G") | starts_with("H"),
               names_to = "Sample",
               values_to = "Abundance") %>%
  mutate(Sample = str_replace(Sample, "Control", "Non-Context")) %>%
  group_by(SuperClass, Sample) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = 'drop')

# Adjust the sample_order for renamed samples
sample_order <- c('A1_Non-Context', 'A1_Context', 'A2_Non-Context', 'A2_Context', 'B1_Non-Context', 'B1_Context', 'B2_Non-Context', 'B2_Context', 'C1_Non-Context', 'C1_Context', 'C2_Non-Context', 'C2_Context', 'D1_Non-Context', 'D1_Context', 'D2_Non-Context', 'D2_Context', 'E1_Non-Context', 'E1_Context', 'E2_Non-Context', 'E2_Context', 'F1_Non-Context', 'F1_Context', 'F2_Non-Context', 'F2_Context', 'G1_Non-Context', 'G1_Context', 'H1_Non-Context', 'H1_Context', 'H25361_Non-Context', 'H25361_Context', 'H25362_Non-Context', 'H25362_Context', 'H25363_Non-Context', 'H25363_Context', 'H25364_Non-Context', 'H25364_Context', 'H25365_Non-Context', 'H25365_Context')
data_sum$Sample <- factor(data_sum$Sample, levels = sample_order)
# Z-score scaling for abundance
data_sum$ZScoreAbundance <- scale(data_sum$Abundance)

# Excluding non-numeric columns from the pivot
metabolite_ttest_results <- data %>%
  select(Metabolite, everything(), -SuperClass, -Class, -SubClass) %>%
  pivot_longer(cols = -Metabolite,
               names_to = "Sample",
               values_to = "Abundance") %>%
  group_by(Metabolite) %>%
  do(tidy = broom::tidy(t.test(.data$Abundance[grepl("_Control$", .data$Sample)],
                               .data$Abundance[grepl("_Context$", .data$Sample)],
                               paired = TRUE))) %>%
  unnest(cols = tidy) %>%
  select(Metabolite, p.value)

# Sorting by the absolute difference of p-values from 0.5 will allow the smallest and largest p-values to be at the top.
# This helps to find the most distinct or similar metabolites.
metabolite_ttest_results <- metabolite_ttest_results %>%
  arrange(abs(p.value - 0.5))

# Display the top 10 metabolites with strongest differences or similarities
head(metabolite_ttest_results, 20)

# Adjust p-values for FDR
metabolite_ttest_results$adjusted_p.value <- p.adjust(metabolite_ttest_results$p.value, method = "BH")

# Display the top metabolites with strongest differences or similarities after FDR correction
head(metabolite_ttest_results, 10)

sorted_data <- metabolite_ttest_results %>%
  arrange(adjusted_p.value)

# Plotting
library(ggplot2)
library(dplyr)

# Assuming sorted_data is already prepared
threshold <- -log10(0.1)  # Define a threshold for significance

# Filtering for significant metabolites (Supp5C)
significant_data <- sorted_data %>%
  filter(-log10(adjusted_p.value) >= threshold)

# Further narrowing less significant metabolites for Supp5D
# Adjust this threshold according to your preference for how much to filter out
less_significant_data_narrowed <- sorted_data %>%
  filter(-log10(adjusted_p.value) < threshold, -log10(adjusted_p.value) > threshold - 2)

# Function to create and save plots with enhanced labeling
create_plot <- function(data, filename) {
  p <- ggplot(data, aes(x = reorder(Metabolite, -adjusted_p.value), y = -log10(adjusted_p.value))) +
    geom_bar(stat = 'identity', aes(fill = -log10(adjusted_p.value))) +
    geom_hline(yintercept = threshold, color = "red", linetype = "dashed") +
    labs(y = "-log10(Adjusted P-Value)", x = "Metabolite") +
    scale_fill_gradientn(colors = c("blue", "yellow", "red"), name = "-log10(Adjusted P-Value)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "bold", size = 7),
          axis.text.y = element_text(face = "bold", size = 14),
          axis.title.x = element_text(size = 14, face = "bold"),
          axis.title.y = element_text(size = 14, face = "bold"),
          legend.position = "bottom")
  
  # Adjusted dimensions for a more box-shaped plot
  ggsave(filename, p, width = 8, height = 6, units = "in", dpi = 300)
}

#Font size 14
# Create and save the plot for significant metabolites (Supp5C)
create_plot(significant_data, "Plot1.png")

# Correcting the extraction of metabolites

# Filter significant metabolites
significant_metabolites <- metabolite_ttest_results %>% 
  filter(adjusted_p.value < 0.05) %>% 
  arrange(adjusted_p.value)

# Identify metabolites with strongest differences
strongest_differences <- head(significant_metabolites, 5)

# Identify metabolites with least differences (i.e., most similarity)
most_similar <- tail(significant_metabolites, 5)

# Print results for verification
print(strongest_differences)
print(most_similar)


# Plotting a heatmap
p <- ggplot(data_sum, aes(x = Sample, y = SuperClass)) +
  geom_tile(aes(fill = ZScoreAbundance), color = "white") +
  scale_fill_gradientn(colors = c("blue", "yellow", "red")) +
  labs(y = "SuperClass", fill = "Z-Score Scaled\nAbundance") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "bold", size = 12, margin = margin(b = 10)), # Rotate x-axis labels to 90 degrees, increase font size, make bold, and add bottom margin
    axis.text.y = element_text(face = "bold", size = 12),  # Increase font size and make y-axis labels bold
    axis.title.y = element_text(face = "bold", size = 14),  # Make y-axis title bold and increase font size
    axis.title.x = element_blank(),  # Remove x-axis title for cleaner look
    axis.ticks.x = element_blank(),  # Remove x-axis ticks for a cleaner look
    plot.margin = unit(c(1, 1, 1, 1), "cm")  # Adjust plot margins to ensure labels fit
  )

# Display the plot
print(p)
ggsave("Plot3.png", plot = p, width = 12, height = 6, units = "in")

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)

# Load the data
data <- read.csv("Biomass.csv")

# Provided sample names
sample_names <- c("A1", "A2", "B1", "B2", "C1", "C2", "D1", "D2",
                  "E1", "E2", "F1", "F2", "G1", "H1", "H25361",
                  "H25362", "H25363", "H25364", "H25365")

# Let's assume samples starting with "H" or containing "2536" are considered "E.coli", others are "L.iners"

# Convert to longer format and extract community and condition
library(ggplot2)
library(dplyr)
library(tidyr)

data_long <- data %>%
  pivot_longer(cols = -Metabolites, names_to = "Sample", values_to = "Flux") %>%
  mutate(
    Community = ifelse(grepl("^H|2536", Sample), "E.coli", "L.iners"),
    Condition = ifelse(str_detect(Sample, "Control"), "Non-Context", "Context") # Rename "Control" to "Non-Context"
  )

# Plot with updated Condition names
p <- ggplot(data_long, aes(x = interaction(Community, Condition, sep = "\n"), y = Flux, fill = Condition)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8), alpha = 0.5) + 
  geom_jitter(width = 0.15, alpha = 0.7, size = 2.5, aes(color = Condition)) + 
  scale_color_manual(values = c("Non-Context" = "#FC4E2A", "Context" = "#E6AB02")) + 
  scale_fill_manual(values = c("Non-Context" = "#FC4E2A", "Context" = "#E6AB02")) + 
  labs(y = "Global Flux", x = "", title = "Comparison of Global Flux between Communities and Conditions") + 
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.title = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.y = element_text(size = 14)
  ) +
  geom_signif(comparisons = list(c("E.coli\nNon-Context", "L.iners\nNon-Context"), c("E.coli\nContext", "L.iners\nContext")),
              map_signif_level = TRUE, textsize = 4, vjust = -0.5,
              # Specify annotations manually
              annotations = c("**", "ns"))  # Adjust annotations as per your data or remove if not applicable

# Print the plot
print(p)
# Assuming 'p' is your ggplot object
ggsave(filename = "Comparative_Flux.png", plot = p, width = 10, height = 8, dpi = 300)
