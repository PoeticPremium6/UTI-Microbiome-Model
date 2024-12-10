library(tidyverse)

# Set working directory to where the .RDS files are stored
setwd("Bacarena/")

# List all .RDS files (corrected to use the uppercase .RDS extension)
files <- list.files(pattern = "\\.RDS$")

# Read the metabolite mapping file with human-readable names
metabolite_mapping <- read_csv("metabolites_mapping.csv")

# Initialize a list to store the average flux data for each organism
avg_flux_data <- list()

# Process each .RDS file
for (file in files) {
  # Read the .RDS file
  sim <- readRDS(file)
  
  # Extract the flux data (assuming it's in sim@mfluxlist as per your previous check)
  flux_data <- sim@mfluxlist
  organism_list <- flux_data[[1]]
  
  # Process each organism's flux data
  for (organism in names(organism_list)) {
    organism_fluxes <- lapply(flux_data, function(dataset) {
      if (!is.null(dataset[[organism]])) {
        ex_reactions <- dataset[[organism]][grep("^EX_", names(dataset[[organism]]))]
        return(ex_reactions)
      }
      return(rep(NA, length(grep("^EX_", names(dataset[[organism]])))))  # Handle missing data
    })
    
    # Calculate average flux for each EX_ reaction for the organism
    avg_flux <- colMeans(do.call(rbind, organism_fluxes), na.rm = TRUE)
    
    # Store the average flux data in the list
    if (!is.null(avg_flux_data[[organism]])) {
      avg_flux_data[[organism]] <- avg_flux_data[[organism]] + avg_flux  # Sum the flux if same organism appears
    } else {
      avg_flux_data[[organism]] <- avg_flux
    }
  }
}

# Convert the avg_flux_data list into a data frame
flux_df <- bind_rows(lapply(names(avg_flux_data), function(organism) {
  tibble(Organism = organism, 
         Metabolite = names(avg_flux_data[[organism]]), 
         Flux = avg_flux_data[[organism]])
}), .id = "File")

# Classify as production (positive flux) or consumption (negative flux)
flux_df <- flux_df %>%
  mutate(Status = case_when(
    Flux > 0 ~ "Production",  # Positive flux indicates production
    Flux < 0 ~ "Consumption",  # Negative flux indicates consumption
    TRUE ~ "None"  # No data or zero flux
  ))

# Calculate the total flux (sum of absolute flux values) for each metabolite across all organisms
total_flux_by_metabolite <- flux_df %>%
  group_by(Metabolite) %>%
  summarize(TotalFlux = sum(abs(Flux), na.rm = TRUE)) %>%
  arrange(desc(TotalFlux))

# Select the top 50 metabolites based on total flux
top_metabolites <- total_flux_by_metabolite$Metabolite[1:50]

# Filter the flux data to include only the top 50 metabolites
flux_df_filtered <- flux_df %>%
  filter(Metabolite %in% top_metabolites)

# Merge the flux data with the human-readable names from the mapping file
flux_df_filtered <- flux_df_filtered %>%
  left_join(metabolite_mapping, by = c("Metabolite" = "Metabolite_ID")) %>%
  mutate(Metabolite = ifelse(is.na(Human_Readable_Name), Metabolite, Human_Readable_Name)) %>%
  select(-Human_Readable_Name)

# Plotting the heatmap with the human-readable names
combined_plot <- ggplot(flux_df_filtered, aes(x = Metabolite, y = Organism, fill = Status)) +
  geom_tile() +
  scale_fill_manual(values = c("Production" = "#4B0082", "Consumption" = "#D3B8E6", "None" = "white")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
  labs(fill = "Flux Status", x = "Metabolite", y = "Organism", title = "")


ggsave("Prod_Cons_Heatmap.png", 
       combined_plot, width = 12, height = 8, dpi = 300)

########################
########################
# Load necessary libraries
library(BacArena)
library(dplyr)
library(igraph)
library(readr)

# List of all sample names
samples <- c("A01", "A02", "B01", "B02", "C01", "C02", "D01", "D02", 
             "E01", "E02", "F01", "F02", "G01", "H01", "H25361", "H25362", 
             "H25363", "H25364", "H25365")

# Initialize an empty list to store all feeding data
all_feeding_data <- list()

# Loop through each sample and extract its feeding data
for (sample in samples) {
  # Load the simulation data for the current sample
  simulation_file <- paste0("Bacarena\\", sample, "_Simulation.RDS")
  A <- readRDS(simulation_file)
  
  # Get a list of the metabolites for the feeding analysis (top 30)
  subs <- names(head(getVarSubs(A), 30))
  
  # Use suppressWarnings() to suppress warnings for samples with no crossfeeding
  feeding_data <- suppressWarnings(findFeeding3(A, time=3, mets=subs, plot=F))
  
  # Check if feeding data is empty or null (i.e., no crossfeeding found)
  if (nrow(feeding_data) == 0) {
    cat("No crossfeeding found for sample", sample, "- Skipping this sample.\n")
    next  # Skip to the next sample in the loop
  }
  
  # Convert the output to a data frame
  feeding_data_matrix <- as.data.frame(feeding_data)
  
  # Add a column to indicate the sample
  feeding_data_matrix$Sample <- sample
  
  # Save the data to a CSV file for the current sample
  output_file <- paste0("\\Bacarena", sample, "_FeedingData.csv")
  write.csv(feeding_data_matrix, output_file, row.names = FALSE)
  
  # Add the data to the list of all feeding data
  all_feeding_data[[sample]] <- feeding_data_matrix
}

# Merge all data frames into one large data frame
merged_feeding_data <- do.call(rbind, all_feeding_data)

# Save the merged data to a CSV file
merged_output_file <- "Bacarena\\Merged_FeedingData.csv"
write.csv(merged_feeding_data, merged_output_file, row.names = FALSE)

cat("All samples processed and merged successfully!")

# Read the metabolite mapping file
metabolite_mapping <- read_csv("metabolites_mapping.csv")

# Check for issues in the mapping file
problems(metabolite_mapping)

# Process the merged feeding data
processed_data <- merged_feeding_data %>%
  # Join with metabolite mapping for human-readable names
  left_join(metabolite_mapping, by = c("met" = "Metabolite_ID")) %>%
  mutate(
    total_flux = abs(prod.flux) + abs(cons.flux)  # Compute total flux as the sum of absolute values
  ) %>%
  group_by(Sample, Human_Readable_Name) %>%  # Group by sample and metabolite
  summarise(
    total_flux_sum = sum(total_flux, na.rm = TRUE)  # Sum total flux for each metabolite per sample
  ) %>%
  ungroup() %>%  # Ungroup for further operations
  arrange(desc(total_flux_sum))  # Sort by highest flux sum

# Save the processed data for visualization
output_file <- "Processed_FeedingData.csv"
write.csv(processed_data, output_file, row.names = FALSE)

cat("Processed data saved successfully!")

# Load ggplot2 for visualization
library(ggplot2)

# Apply log normalization to the total_flux_sum
processed_data <- processed_data %>%
  mutate(log_total_flux = log1p(total_flux_sum))  # log1p adds 1 to avoid log(0)

custom_palette <- c(
  "#D55E00", "#E69F00", "#F0E442", "#009E73", "#56B4E9",
  "#0072B2", "#CC79A7", "#999999", "#FF6666", "#77AC30",
  "#EDB120", "#7E2F8E", "#4DBEEE", "#A2142F", "#6666FF",
  "#FF00FF", "#00FFFF", "#CCFF00", "#FF9999"
)

# Create the stacked bar plot with black borders between the bars
combined_plot <- ggplot(processed_data, aes(x = Human_Readable_Name, y = log_total_flux, fill = Sample)) +
  geom_bar(stat = "identity", color = "black") +  # Add black border around each bar
  scale_fill_manual(values = custom_palette) +  # Apply the custom color palette
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14, face = "bold"),  # Increase size and bold
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 15, face = "bold"),  # Increase size for x-axis title
    axis.title.y = element_text(size = 15, face = "bold"),  # Increase size for y-axis title
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  labs(
    fill = "Sample",
    x = "Metabolite",
    y = "Log Normalized Rate of Metabolite Exchange",
    title = ""
  )

# Save the plot as a high-resolution PNG
ggsave(
  filename = "/Log_Normalized_Flux_Stacked_Bar.png",
  plot = combined_plot,
  width = 16,
  height = 10,
  dpi = 300
)
