library(BacArena)
library(dplyr)
library(igraph)

# Define a function to process a single file
process_simulation_file <- function(file_path, output_dir) {
  cat("Processing file:", file_path, "\n")
  sim <- readRDS(file_path)
  
  # Ensure metabolites are correctly identified
  if (is.null(sim@mediac) || length(names(sim@mediac)) == 0) {
    cat("No metabolites found in sim@mediac for file:", file_path, "\n")
    return()
  }
  
  mets <- names(sim@mediac)
  
  # Apply the analysis
  g <- findFeeding3(sim, time = 4, mets = sim@mediac, cutoff = 1e-12)
  
  if (is.null(g) || length(g) == 0) {
    cat("No cross-feeding interactions found for file:", file_path, "\n")
    return()
  }
  
  # Process and format the output
  combined_df <- do.call(rbind, g)
  
  formatted_df <- combined_df %>%
    mutate(
      sample = gsub(".*/|\\.rds$", "", file_path),  # Extracting sample name from file path
      producer = gsub("\\.mat$", "", prod),  # Removing the .mat extension from producer names
      consumer = gsub("\\.mat$", "", cons),  # Removing the .mat extension from consumer names
      metabolite = gsub("EX_", "", met),  # Simplifying metabolite ID
      metabolite = gsub("_e0$", "", metabolite),  # Further simplifying metabolite ID
      interaction_type = case_when(
        prod.flux > 0 & cons.flux < 0 ~ "Produced",  # Producer if positive prod.flux and negative cons.flux
        prod.flux < 0 & cons.flux > 0 ~ "Consumed",  # Consumer if negative prod.flux and positive cons.flux
        TRUE ~ "Neutral"  # Neutral for other cases
      )
    ) %>%
    select(sample, producer, consumer, metabolite, interaction_type, prod.flux, cons.flux) %>%  # Including prod.flux and cons.flux
    filter(interaction_type != "Neutral")  # Optional: filter out neutral interactions
  
  
  # Define the output file name based on the input file's name
  output_file_path <- paste0(output_dir, "/", gsub(".*/|\\.rds$", "", file_path), "_cross_feeding_interactions.csv")
  
  # Save the results
  write.csv(formatted_df, output_file_path, row.names = FALSE, quote = FALSE)
  cat("Saved cross-feeding interactions to:", output_file_path, "\n")
}

# Directory paths
input_dir <- "/UTI/Crossfeed"
output_dir <- "/UTI/Figures/Crossfeeding"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Get a list of all .RDS files
sim_files <- list.files(input_dir, full.names = TRUE, pattern = "\\.rds$")

# Process each .RDS file
lapply(sim_files, process_simulation_file, output_dir)

#Now that we have the data, let's merge them together
library(dplyr)
library(readr)

# Create the output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Set the working directory to where your files are located
setwd(input_dir)

# List all CSV files in the directory
files <- list.files(pattern = "_cross_feeding_interactions\\.csv$")

# Function to extract the base sample name (e.g., A1_Context, A1_Control) from the filename
extract_sample_name <- function(filename) {
  sub("(_\\d+)?_cross_feeding_interactions\\.csv$", "", filename)
}

# Group files by their base sample name
grouped_files <- split(files, sapply(files, extract_sample_name))

# Loop through each group and merge files, including the sample name in the merged data
lapply(names(grouped_files), function(sample_name) {
  # For each file in the group, read the file, add a 'SampleName' column, and combine
  combined_df <- lapply(grouped_files[[sample_name]], function(file) {
    df <- read_csv(file)
    # Add a new column with the base sample name
    mutate(df, SampleName = sample_name)
  }) %>% bind_rows()
  
  # Write the combined data to a new master file in the output directory
  master_file_path <- paste0(output_dir, "/", sample_name, "_cross_feeding_interactions_master.csv")
  write_csv(combined_df, master_file_path)
  
  cat("Merged files for", sample_name, "into", master_file_path, "\n")
})
#Let's clean up the metabolite naming
library(dplyr)
library(readr)

# Define the directory containing the merged interaction data files
data_dir <- "/Users/josspa/UTI/Figures/Crossfeeding/Merged"

# Define the path to the mapping file
mapping_file_path <- paste0(data_dir, "/unique_metabolites_clean.csv")

# Read the mapping file into a dataframe
metabolite_mapping <- read_csv(mapping_file_path)

# List all interaction data CSV files in the directory
interaction_files <- list.files(path = data_dir, pattern = "A1_Context_cross_feeding_interactions_master\\.csv$", full.names = TRUE)

# Function to process and save each interaction data file with human-readable metabolite names
process_and_save_files <- function(file_path, mapping) {
  # Read interaction data
  interaction_data <- read_csv(file_path)
  
  # Replace metabolite IDs with human-readable names using a left join
  updated_data <- interaction_data %>%
    left_join(mapping, by = c("metabolite" = "Metabolite_ID")) %>%
    mutate(metabolite = coalesce(Human_Readable_Name, metabolite)) %>%
    select(-Human_Readable_Name)
  
  # Save the updated interaction data back to the same file
  write_csv(updated_data, file_path)
  
  message("Processed and saved: ", basename(file_path))
}

# Apply the processing function to each file
lapply(interaction_files, process_and_save_files, metabolite_mapping)

#Now we can plot everything!
library(tidyverse)

# Read and combine all files into one dataframe
combined_data <- list.files(path = input_dir, pattern = "\\.csv$", full.names = TRUE) %>%
  map_dfr(read_csv) %>%
  mutate(sample = str_remove(sample, "_\\d+$"))  # Group similar samples

library(tidyverse)

# Filter for 'Context' samples
context_samples_summary <- combined_data %>%
  filter(grepl("Context$", SampleName)) %>%
  count(SampleName, metabolite, name = "frequency")

# Identify top 20 metabolites based on total frequency
top_20_metabolites <- context_samples_summary %>%
  group_by(metabolite) %>%
  summarise(total_frequency = sum(frequency), .groups = 'drop') %>%
  arrange(desc(total_frequency)) %>%
  slice_head(n = 20) %>%
  pull(metabolite)

# Filter the dataset for these top 20 metabolites
top_20_context_data <- context_samples_summary %>%
  filter(metabolite %in% top_20_metabolites)

# Prepare a color palette
n_samples <- length(unique(top_20_context_data$SampleName))
palette <- brewer.pal(min(n_samples, 8), "Dark2")
if (n_samples > 8) {
  palette <- colorRampPalette(brewer.pal(8, "Dark2"))(n_samples)
}

# Create the plot for top 20 metabolite exchanges
ggplot(top_20_context_data, aes(x = reorder(metabolite, -frequency), y = frequency, fill = SampleName)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +  # Flipping coordinates for better visualization of metabolites
  scale_fill_manual(values = palette) +  # Custom palette for SampleName
  labs(title = "", 
       x = "Rate of Metabolite Exchange", 
       y = "Metabolite", 
       fill = "Sample Name") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),  # Centered, bold, and larger title
    legend.title = element_text(face = "bold", size = 14),  # Bold and larger legend title
    legend.text = element_text(face = "bold", size = 12),  # Larger legend text
    legend.position = "bottom",  # Legend at the bottom
    axis.title.x = element_text(face = "bold", size = 14),  # Bold and larger X-axis title
    axis.title.y = element_text(face = "bold", size = 14),  # Bold and larger Y-axis title
    axis.text.x = element_text(face = "bold", size = 12),  # Bold and larger X-axis text
    axis.text.y = element_text(face = "bold", size = 12)  # Bold and larger Y-axis text
  )

# Save the plot
ggsave("Top20_Metabolite_Exchanges_Context_Samples.png", width = 12, height = 10, units = "in")
library(ggplot2)
library(dplyr)
library(scales) # For additional color scales if needed
library(RColorBrewer) # For color palettes

# Assuming combined_data exists and is correctly formatted

# Step 1: Filter for 'Context' samples and count frequencies
context_samples_summary <- combined_data %>%
  filter(grepl("Context$", SampleName)) %>%
  count(SampleName, metabolite, name = "frequency")

# Step 2: Identify top 20 metabolites based on total frequency
top_20_metabolites <- context_samples_summary %>%
  group_by(metabolite) %>%
  summarise(total_frequency = sum(frequency), .groups = 'drop') %>%
  arrange(desc(total_frequency)) %>%
  slice_head(n = 20) %>%
  pull(metabolite)

# Step 3: Filter the dataset for these top 20 metabolites
top_20_context_data <- context_samples_summary %>%
  filter(metabolite %in% top_20_metabolites)

# Ensure steps 1 to 3 have been executed without error before proceeding.

# Step 4: Prepare a color palette
n_samples <- length(unique(top_20_context_data$SampleName))
palette <- brewer.pal(min(n_samples, 8), "Dark2")
if (n_samples > 8) {
  palette <- colorRampPalette(brewer.pal(8, "Dark2"))(n_samples)
}

# Step 5: Create the plot for top 20 metabolite exchanges, with improvements
plot <- ggplot(top_20_context_data, aes(x = reorder(metabolite, -frequency), y = frequency, fill = SampleName)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
  scale_fill_manual(values = palette, name = "Sample Name") +
  labs(title = "Top 20 Metabolite Exchanges in Context Samples", x = "Rate of Metabolite Exchange", y = "Metabolite") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(face = "bold", size = 12),
        legend.position = "right",
        axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 12),
        axis.text.y = element_text(face = "bold", size = 12))

# Display the plot
print(plot)
ggsave("Top20_Metabolite_Exchanges_Context_Samples.png", width = 12, height = 10, units = "in")

# Prepare detailed information for a CSV file
# Assuming you want detailed interaction info for the top 20 metabolites across all 'Context' samples

# Filter original combined_data for top 20 metabolites (not just summary counts)
detailed_context_data_top_20 <- combined_data %>%
  filter(metabolite %in% top_20_metabolites, grepl("Context$", SampleName))

# Save this detailed interaction data to a CSV file
write_csv(top_20_context_data, "Top20_Metabolites_Context_Interactions.csv")

library("BacArena")
library("R.matlab")
library("parallel")

#Metabolites of interest
#Alanine (EX_cpd00035_e0)
#ornithine (EX_cpd00064_e0)
#Acetic Acid (EX_cpd00029_e0)
#L-Serine (EX_cpd00054_e0)
#Glycerol (EX_cpd00100_e0)
#Fe2_ (EX_cpd10515_e0)
#Fe3_ (EX_cpd10516_e0)
#Sorbitol (EX_cpd00588_e0)
#Citrate (EX_cpd00137_e0)
#Succinate (EX_cpd00036_e0)
#Butyrate (EX_cpd00211_e0)
#glycine (EX_cpd00033_e0)
#Propionate (EX_cpd00141_e0)
#Mannitol (EX_cpd00314_e0)
#Xylitol (EX_cpd00306_e0)
sim <- readRDS(file = "A1_Context.rds")
#plotSpecActivity(simlist = sim, useNames = T)

#Alanine (EX_cpd00035_e0)
#ornithine (EX_cpd00064_e0)
#H1 is an interesting case
#E.coli is providing ornithine to numerous microbes and E. coli is also receiving this metabolite
#E.coli is also providing sorbitol to 3 other members of the consortium
#H25364 is another interesting case
#3 known gut-associated microbes are receiving both Ornithinine and iron
#More so, E.coli is the core producer of mannitol, an alcohol sugar into the community


g <- findFeeding3(sim, time = 4, mets = c("EX_cpd00064_e0", "EX_cpd10515_e0"), cutoff = 1e-12)

#Mannitol (EX_cpd00314_e0)
#Glycerol (EX_cpd00100_e0)
#Xylitol (EX_cpd00306_e0)
#Sorbitol (EX_cpd00588_e0)
g <- findFeeding3(sim, time = 4, mets = c("EX_cpd00100_e0", 
                                          "EX_cpd00588_e0", 
                                          "EX_cpd00314_e0", 
                                          "EX_cpd00306_e0"), cutoff = 1e-12)

png(filename = "UTI/Figures/crossfeeding", width = 800, height = 800)
findFeeding3(sim, time = 4, mets = c("EX_cpd00035_e0", "EX_cpd00064_e0", "EX_cpd00100_e0"), cutoff = 1e-3)
dev.off()
