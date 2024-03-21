#Producing/Consuming Matrix:
library(tidyverse)
library(ggplot2)

# Set working directory to where the .RDS files are stored
setwd("rUTI/Crossfeed/")

# Read the metabolites mapping file
metabolites_mapping <- read_csv("metabolites_mapping.csv")

# List all .RDS files
files <- list.files(pattern = "\\.rds$")

# Initialize a list to store the average flux data for each organism
avg_flux_data <- list()

# Process each .RDS file
for (file in files) {
  sim <- readRDS(file)
  flux_data <- sim@mfluxlist
  organism_list <- flux_data[[1]]
  
  # Process each organism's flux data
  for (organism in names(organism_list)) {
    organism_fluxes <- lapply(flux_data, function(dataset) {
      if (!is.null(dataset[[organism]])) {
        ex_reactions <- dataset[[organism]][grep("^EX_", names(dataset[[organism]]))]
        return(ex_reactions)
      }
      return(rep(NA, length(grep("^EX_", names(dataset[[organism]])))))
    })
    
    # Calculate average flux for each EX_ reaction for the organism
    avg_flux <- colMeans(do.call(rbind, organism_fluxes), na.rm = TRUE)
    avg_flux_data[[organism]] <- avg_flux
  }
}

# Clean organism names, remove '.mat', '.mat_1', '_control', and correct spelling
cleaned_avg_flux_data <- list()
for (organism in names(avg_flux_data)) {
  cleaned_name <- gsub("(\\.mat(_[0-9]+)?$)|(_control$)", "", organism)
  cleaned_name <- gsub("Eschericia_coli", "Escherichia_coli", cleaned_name)
  
  if (!is.null(cleaned_avg_flux_data[[cleaned_name]])) {
    cleaned_avg_flux_data[[cleaned_name]] <- (cleaned_avg_flux_data[[cleaned_name]] + avg_flux_data[[organism]]) / 2
  } else {
    cleaned_avg_flux_data[[cleaned_name]] <- avg_flux_data[[organism]]
  }
}

# Prepare combined_flux_df with corrected organism names
combined_flux_df <- map_df(cleaned_avg_flux_data, ~tibble(Metabolite = names(.x), Flux = .x), .id = "Organism")

# Correct organism names in the dataframe
combined_flux_df$Organism <- gsub("(\\.mat(_[0-9]+)?$)|(_control$)", "", combined_flux_df$Organism)
combined_flux_df$Organism <- gsub("Eschericia_coli", "Escherichia_coli", combined_flux_df$Organism)

# Join with metabolites mapping
combined_flux_df <- combined_flux_df %>%
  left_join(metabolites_mapping, by = c("Metabolite" = "Metabolite_ID")) %>%
  mutate(Metabolite = ifelse(is.na(Human_Readable_Name), Metabolite, Human_Readable_Name),
         Status = case_when(
           Flux > 0 ~ "Produced",
           Flux < 0 ~ "Consumed"
         )) %>%
  filter(Status %in% c("Produced", "Consumed")) # Remove rows where Status would be "Not Applicable"

# For better x-axis readability, let's filter for a subset of metabolites with the highest variance in Flux
top_metabolites <- combined_flux_df %>%
  group_by(Metabolite) %>%
  summarize(Variance = var(Flux, na.rm = TRUE)) %>%
  top_n(50, Variance) %>%
  pull(Metabolite)

combined_flux_df_filtered <- combined_flux_df %>%
  filter(Metabolite %in% top_metabolites)

# Plotting with flipped y-axis and filtered metabolites
ggplot(combined_flux_df_filtered, aes(x = Metabolite, y = fct_rev(Organism), fill = Status)) +
  geom_tile() +
  scale_fill_manual(values = c("Consumed" = "blue", "Produced" = "green")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, face = "bold"), 
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
  labs(fill = "Status", x = "Metabolite", y = "Organism", title = "Metabolite Consumption and Production")

# Save the plot in landscape orientation
ggsave(paste0(output_dir, "Filtered_Metabolite_Status_Heatmap_150_Landscape.png"), plot = last_plot(), width = 25, height = 8, dpi = 300)


# Ensure output directory exists and save the plot
output_dir <- "/Users/josspa/UTI/Figures/Figure6/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
ggsave(paste0(output_dir, "Filtered_Metabolite_Status_Heatmap_50.png"), plot = last_plot(), width = 20, height = 8, dpi = 300)

library(ggplot2)
library(dplyr)
library(patchwork)

# Subset for Produced metabolites
produced_df <- combined_flux_df_filtered %>%
  filter(Status == "Produced")

# Subset for Consumed metabolites
consumed_df <- combined_flux_df_filtered %>%
  filter(Status == "Consumed")

# Plot for Produced metabolites
plot_produced <- ggplot(produced_df, aes(x = Metabolite, y = fct_rev(Organism), fill = Status)) +
  geom_tile() +
  scale_fill_manual(values = c("Produced" = "green")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 20, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "none"  # Remove legend if it's not necessary
  ) +
  labs(x = "", y = "Organism", title = "Metabolites Produced")

# Improve the plot for Consumed metabolites
plot_consumed <- ggplot(consumed_df, aes(x = Metabolite, y = fct_rev(Organism), fill = Status)) +
  geom_tile() +
  scale_fill_manual(values = c("Consumed" = "blue")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 20, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "none"  # Remove legend if it's not necessary
  ) +
  labs(x = "", y = "Organism", title = "Metabolites Consumed")

# Print the improved plots to check the adjustments
print(plot_produced)
print(plot_consumed)

# Combine the plots
combined_plot <- plot_produced / plot_consumed

# Display the combined plot
print(combined_plot)

# Save the combined plot
ggsave(paste0(output_dir, "Metabolite_Status_Split_Heatmap_50_new.png"), plot = combined_plot, width = 25, height = 20, dpi = 300)

#Now let's subset the previous part for specific organisms
# Define the selected organisms
selected_species <- c("Escherichia_coli", "Lactobacillus_iners", "Shigella_sonnei", "Salmonella_enterica", "Levilactobacillus paucivorans")

# Filter the combined_flux_df_filtered dataframe for only the selected organisms
combined_flux_df_subset <- combined_flux_df_filtered %>%
  filter(Organism %in% selected_species)


library(ggplot2)
library(dplyr)

# Assuming combined_flux_df_filtered and selected_species are already defined and available

# Filter the combined_flux_df_subset dataframe for only the selected organisms
combined_flux_df_subset <- combined_flux_df_filtered %>%
  filter(Organism %in% selected_species)

# Subset for Produced metabolites
produced_df_subset <- combined_flux_df_subset %>%
  filter(Status == "Produced")

# Subset for Consumed metabolites
consumed_df_subset <- combined_flux_df_subset %>%
  filter(Status == "Consumed")


library(forcats) # For fct_rev

# Assuming combined_flux_df_filtered is already defined and available
# Assuming selected_species is defined and combined_flux_df_subset is filtered

# Combine produced and consumed subsets for plotting
combined_df_for_plotting <- rbind(produced_df_subset, consumed_df_subset)

# Plot both Produced and Consumed metabolites in a single plot
combined_plot <- ggplot(combined_df_for_plotting, aes(x = Metabolite, y = fct_rev(Organism), fill = Status)) +
  geom_tile() +
  scale_fill_manual(values = c("Produced" = "green", "Consumed" = "blue")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14, face = "bold"),
        axis.text.y = element_text(size = 18, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
  labs(fill = "Metabolic Status", x = "Metabolite", y = "Organism", title = "Metabolite Production and Consumption by Selected Microbes")

# Display the combined plot
print(combined_plot)

# Save the plots
output_dir <- "/Figure7"
# Ensure the output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save the produced plot
ggsave(paste0(output_dir, "Selected_Species_Heatmap.png"), plot = combined_plot, width = 12, height = 12, dpi = 300)
