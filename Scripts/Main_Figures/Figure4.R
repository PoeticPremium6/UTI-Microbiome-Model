###FIGURE 4####
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

# Plotting 
p <- ggplot(data_long, aes(x = Metabolite, y = variable)) + 
  geom_tile(aes(fill = factor(value)), color = "white") + 
  scale_fill_viridis_d(name = "Presence", begin = 0.2, end = 0.8,
                       labels = c("0" = "Absent", "1" = "Present")) + 
  labs(x = "Metabolite", y = "Sample") + 
  theme_minimal(base_size = 14) + 
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"), 
    axis.title = element_text(size = 16, face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 12),
    plot.margin = margin(1, 1, 1, 1, "cm")
  ) +
  coord_fixed(ratio = 1) +
  geom_hline(yintercept = seq(2.5, length(unique(data_long$variable)), by = 2), color = "black", size = 1.5)  # Adjusted line placement

# Display the plot
print(p)

# Save the plot 
ggsave("metabolite_heatmap_differences_filtered.png", p, width = 16, height = 12, units = "in")


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
