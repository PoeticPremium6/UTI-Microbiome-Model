# Load necessary libraries
library(dplyr)  # For data manipulation
library(tidyr)  # For data reshaping
library(ggplot2) # For plotting
library(tibble) # For tidy data manipulation

# Define the input file path
input_file <- "Combined_Reaction_Flux_Summary.csv"

# Step 1: Read the CSV file into a dataframe
data <- read.csv(input_file, stringsAsFactors = FALSE)

# Define human-readable names mapping (example provided by user)
reaction_name_map <- tibble(
  Reaction = c(
    "HXANtex", "HYXNtpp", "7OCHOLATEtex", "7ocholate[u]tr",
    "TDG3PACT", "CADVtpp", 
    "5DGLCNR", "IDON_Ltex", 
    "IDONt2rpp", "idon_L[u]tr", 
    "FDNADOX_Hpp", "CHOLATEtex", 
    "cholate[u]tr", "G3PAT140", 
    "GLCP", "GLCS1", 
    "15DAPtpp", "KARA1",
    "GLCNtex", "glcn[u]tr"
  ),
  HumanReadableName = c(
    "Hypoxanthine Exchange", "Hypoxanthine Transport", "Cholanate Exchange", "Cholanate Transport",
    "Myristoyl O-acyltransferase", "Lysine/Cadaverine antiporter", 
    "5-dehydro-D-gluconate reductase", "L-Idonate exchange", 
    "L-idonate transport via proton symport", "L-Idonate Transport",
    "Ferredoxin", "Cholate Transport via Bicarbonate Antiport",
    "Cholate Transport", "Glycerol-3-phosphate acyltransferase",
    "Glycogen phosphorylase", "Glycogen synthase",
    "1,5-Diaminopentane (cadaverine) export", "Ketol-acid reductoisomerase",
    "D-gluconate transport via proton symport","D-gluconate Transport"
  )
)

#######################################
# Positive Flux Analysis (Production)
# Step 2: Filter for positive values
positive_data <- data %>%
  mutate(across(where(is.numeric), ~ ifelse(. > 0, ., NA)))  # Retain positive values, replace others with NA

# Step 3: Reshape and summarize positive flux data
long_positive_data <- positive_data %>%
  pivot_longer(-Reaction, names_to = "Sample", values_to = "Value") %>%
  filter(!is.na(Value))  # Retain only rows with positive values

# Perform a one-sample t-test for each reaction (positive flux)
stat_results_positive <- long_positive_data %>%
  group_by(Reaction) %>%
  summarise(
    mean_flux = mean(Value, na.rm = TRUE),  # Mean positive flux
    sd_flux = sd(Value, na.rm = TRUE),     # Standard deviation
    n = n(),                               # Sample size
    p_value = ifelse(
      n > 1 & sd_flux > 0,                 # Only perform t-test if more than 1 value and variability exists
      t.test(Value, mu = 0)$p.value,       # One-sample t-test (against null hypothesis of zero flux)
      NA                                  # Assign NA if test cannot be performed
    ),
    .groups = "drop"
  )

# Filter reactions based on p-value and positive mean flux thresholds
significant_positive_reactions <- stat_results_positive %>%
  filter(p_value < 0.05 & mean_flux > 750)  # Adjust threshold as needed (positive mean flux and p-value < 0.05)

# Select the top 10 positive reactions, first by mean flux, then by p-value (highest confidence)
top_positive_reactions <- significant_positive_reactions %>%
  arrange(desc(mean_flux), p_value) %>%
  slice_head(n = 10)

# Step 4: Prepare the positive matrix for plotting (with HumanReadableName)
positive_matrix <- positive_data %>%
  pivot_longer(-Reaction, names_to = "Sample", values_to = "Value") %>%  # Reshape to long format
  filter(Reaction %in% top_positive_reactions$Reaction)  # Retain top positive reactions

# Join with the human-readable names for positive reactions
positive_matrix <- positive_matrix %>%
  left_join(reaction_name_map, by = "Reaction")

positive_matrix$flux_type <- "Production"  # Label as Production

#######################################
# Negative Flux Analysis (Consumption)
# Step 2: Filter for negative values
negative_data <- data %>%
  mutate(across(where(is.numeric), ~ ifelse(. < 0, ., NA)))  # Retain negative values, replace others with NA

# Step 3: Reshape and summarize negative flux data
long_negative_data <- negative_data %>%
  pivot_longer(-Reaction, names_to = "Sample", values_to = "Value") %>%
  filter(!is.na(Value))  # Retain only rows with negative values

# Perform a one-sample t-test for each reaction (negative flux)
stat_results_negative <- long_negative_data %>%
  group_by(Reaction) %>%
  summarise(
    mean_flux = mean(Value, na.rm = TRUE),  # Mean negative flux
    sd_flux = sd(Value, na.rm = TRUE),     # Standard deviation
    n = n(),                               # Sample size
    p_value = ifelse(
      n > 1 & sd_flux > 0,                 # Only perform t-test if more than 1 value and variability exists
      t.test(Value, mu = 0)$p.value,       # One-sample t-test (against null hypothesis of zero flux)
      NA                                  # Assign NA if test cannot be performed
    ),
    .groups = "drop"
  )

# Filter reactions based on p-value and negative mean flux thresholds
significant_negative_reactions <- stat_results_negative %>%
  filter(p_value < 0.05 & mean_flux < -750)  # Adjust threshold as needed (negative mean flux and p-value < 0.05)

# Select the top 10 negative reactions, first by mean flux, then by p-value (highest confidence)
top_negative_reactions <- significant_negative_reactions %>%
  arrange(mean_flux, p_value) %>%
  slice_head(n = 10)

# Step 4: Prepare the negative matrix for plotting (with HumanReadableName)
negative_matrix <- negative_data %>%
  pivot_longer(-Reaction, names_to = "Sample", values_to = "Value") %>%  # Reshape to long format
  filter(Reaction %in% top_negative_reactions$Reaction)  # Retain top negative reactions

# Join with the human-readable names for negative reactions
negative_matrix <- negative_matrix %>%
  left_join(reaction_name_map, by = "Reaction")

negative_matrix$flux_type <- "Consumption"  # Label as Consumption

#######################################
# Combine Positive and Negative Data for Plotting
combined_data <- bind_rows(
  positive_matrix,
  negative_matrix
)

# Step 5: Plot the top 10 positive and top 10 negative fluxes
# Plot with improved font size and contrasting purple color palette
ggplot(combined_data, aes(x = reorder(HumanReadableName, Value), y = Value, fill = flux_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +  # Flip the coordinates for better readability
  labs(
    title = "",
    x = "Reactions",
    y = "Mean Flux Value",
    fill = "Flux Type"
  ) +
  scale_fill_manual(values = c("Production" = "#4B0082", "Consumption" = "#D3B8E6")) +  # Dark purple and pale purple
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 14),  # Increase font size for x axis
    axis.text.y = element_text(face = "bold", size = 14),  # Increase font size for y axis
    axis.title.x = element_text(face = "bold", size = 16),  # Increase font size for x axis title
    axis.title.y = element_text(face = "bold", size = 16),  # Increase font size for y axis title
    legend.text = element_text(face = "bold", size = 14),  # Increase font size for legend text
    legend.title = element_text(face = "bold", size = 16)  # Increase font size for legend title
  )


# Define the output file path
# Define the output file path
output_file <- "Reaction_Flux_Plot.png"

# Save the plot with updated color palette and larger font size
ggsave(output_file, 
       plot = ggplot(combined_data, aes(x = reorder(HumanReadableName, Value), y = Value, fill = flux_type)) +
         geom_bar(stat = "identity", position = "dodge") +
         coord_flip() +  # Flip the coordinates for better readability
         labs(
           title = "",
           x = "Reactions",
           y = "Mean Flux Value",
           fill = "Flux Type"
         ) +
         scale_fill_manual(values = c("Production" = "#4B0082", "Consumption" = "#D3B8E6")) +  # Dark purple and pale purple
         theme_minimal() +
         theme(
           axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 14),  # Increase font size for x axis
           axis.text.y = element_text(face = "bold", size = 14),  # Increase font size for y axis
           axis.title.x = element_text(face = "bold", size = 16),  # Increase font size for x axis title
           axis.title.y = element_text(face = "bold", size = 16),  # Increase font size for y axis title
           legend.text = element_text(face = "bold", size = 14),  # Increase font size for legend text
           legend.title = element_text(face = "bold", size = 16)  # Increase font size for legend title
         ),
       width = 10, height = 7)  # Adjust size if necessary


# Print the top 10 positive reactions
cat("\nTop 10 Positive (Production) Reactions:\n")
print(top_positive_reactions)

# Print the top 10 negative reactions
cat("\nTop 10 Negative (Consumption) Reactions:\n")
print(top_negative_reactions)

# Print a message to confirm plot completion
cat("\nPlot generated for top 10 positive (Production) and top 10 negative (Consumption) flux reactions.\n")

#################################
#################################
#################################
#################################
#################################
#################################

# Load necessary libraries
library(dplyr)  # For data manipulation
library(tidyr)  # For data reshaping
library(ggplot2) # For plotting
library(tibble) # For tidy data manipulation

# Define the input file path
input_file <- "Combined_Metabolite_Uptake_Production_Summary.csv"

# Step 1: Read the CSV file into a dataframe
data <- read.csv(input_file, stringsAsFactors = FALSE)

# Create a mapping for the top metabolites to human-readable names
metabolite_name_map <- tibble(
  Metabolite = c(
    "hxan[p]", "h[c]", "7ocholate[u]", "idon_L[p]", "idon_L[u]", "3hbcoa_R[c]", 
    "dms[u]", "isobut[p]", "isobut[u]", "nadh[c]", 
    "co2[fe]", "nh4[fe]", "h[fe]", "etoh[fe]", "lac_L[fe]", "dad_2[fe]", 
    "glcur[fe]", "ala_D[fe]", "gsn[fe]"
  ),
  HumanReadableName = c(
    "Hypoxanthine (periplasm)", "Hydrogen ion (cytoplasm)", "7-Oxocholate (uptake)", 
    "L-Idonate (periplasm)", "L-Idonate (uptake)", "3-Hydroxybutyryl-CoA (cytoplasm)",
    "Dimethyl sulfide (uptake)", "Isobutyrate (periplasm)", "Isobutyrate (uptake)", 
    "NADH (cytoplasm)", "Carbon dioxide (fermentation)", "Ammonium (fermentation)", 
    "Hydrogen ion (fermentation)", "Ethanol (fermentation)", "L-Lactate (fermentation)", 
    "Deoxyadenosine diphosphate (fermentation)", "Glucuronate (fermentation)", 
    "D-Alanine (fermentation)", "Guanosine (fermentation)"
  )
)

#######################################
# Positive Flux Analysis (Production)
# Step 2: Filter for positive values
positive_data <- data %>%
  mutate(across(where(is.numeric), ~ ifelse(. > 0, ., NA)))

# Step 3: Reshape and summarize positive flux data
long_positive_data <- positive_data %>%
  pivot_longer(-Metabolite, names_to = "Sample", values_to = "Value") %>%
  filter(!is.na(Value))

stat_results_positive <- long_positive_data %>%
  group_by(Metabolite) %>%
  summarise(
    mean_flux = mean(Value, na.rm = TRUE),
    sd_flux = sd(Value, na.rm = TRUE),
    n = n(),
    p_value = ifelse(
      n > 1 & sd_flux > 0,
      t.test(Value, mu = 0)$p.value,
      NA
    ),
    .groups = "drop"
  )

significant_positive_reactions <- stat_results_positive %>%
  filter(p_value < 0.05 & mean_flux > 750)

top_positive_reactions <- significant_positive_reactions %>%
  arrange(desc(mean_flux), p_value) %>%
  slice_head(n = 10)

positive_matrix <- positive_data %>%
  pivot_longer(-Metabolite, names_to = "Sample", values_to = "Value") %>%
  filter(Metabolite %in% top_positive_reactions$Metabolite) %>%
  left_join(metabolite_name_map, by = "Metabolite")

positive_matrix$flux_type <- "Production"

#######################################
# Negative Flux Analysis (Consumption)
negative_data <- data %>%
  mutate(across(where(is.numeric), ~ ifelse(. < 0, ., NA)))

long_negative_data <- negative_data %>%
  pivot_longer(-Metabolite, names_to = "Sample", values_to = "Value") %>%
  filter(!is.na(Value))

stat_results_negative <- long_negative_data %>%
  group_by(Metabolite) %>%
  summarise(
    mean_flux = mean(Value, na.rm = TRUE),
    sd_flux = sd(Value, na.rm = TRUE),
    n = n(),
    p_value = ifelse(
      n > 1 & sd_flux > 0,
      t.test(Value, mu = 0)$p.value,
      NA
    ),
    .groups = "drop"
  )

significant_negative_reactions <- stat_results_negative %>%
  filter(p_value < 0.05 & mean_flux < -750)

top_negative_reactions <- significant_negative_reactions %>%
  arrange(mean_flux, p_value) %>%
  slice_head(n = 10)

negative_matrix <- negative_data %>%
  pivot_longer(-Metabolite, names_to = "Sample", values_to = "Value") %>%
  filter(Metabolite %in% top_negative_reactions$Metabolite) %>%
  left_join(metabolite_name_map, by = "Metabolite")

negative_matrix$flux_type <- "Consumption"

#######################################
# Combine and Plot
# Join top reactions with human-readable names
top_positive_reactions <- top_positive_reactions %>%
  left_join(metabolite_name_map, by = "Metabolite")

top_negative_reactions <- top_negative_reactions %>%
  left_join(metabolite_name_map, by = "Metabolite")

# Check for unmatched metabolites
unmatched_positives <- setdiff(top_positive_reactions$Metabolite, metabolite_name_map$Metabolite)
unmatched_negatives <- setdiff(top_negative_reactions$Metabolite, metabolite_name_map$Metabolite)

if (length(unmatched_positives) > 0 || length(unmatched_negatives) > 0) {
  cat("\nWarning: The following metabolites were not found in the mapping table:\n")
  cat("Unmatched positives: ", unmatched_positives, "\n")
  cat("Unmatched negatives: ", unmatched_negatives, "\n")
}

# Prepare combined data for plotting
positive_matrix <- positive_matrix %>%
  left_join(metabolite_name_map, by = "Metabolite")

negative_matrix <- negative_matrix %>%
  left_join(metabolite_name_map, by = "Metabolite")

combined_data <- bind_rows(
  positive_matrix %>% filter(!is.na(HumanReadableName)),
  negative_matrix %>% filter(!is.na(HumanReadableName))
)

# Plot the data with enhanced readability and purple color palette
ggplot(combined_data, aes(x = reorder(HumanReadableName, Value), y = Value, fill = flux_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(
    title = "",
    x = "Metabolites",
    y = "Mean Flux Value",
    fill = "Flux Type"
  ) +
  scale_fill_manual(values = c("Production" = "#4B0082", "Consumption" = "#D3B8E6")) +  # Dark purple and pale purple
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 14),  # Increase font size for x axis
    axis.text.y = element_text(face = "bold", size = 14),  # Increase font size for y axis
    axis.title.x = element_text(face = "bold", size = 16),  # Increase font size for x axis title
    axis.title.y = element_text(face = "bold", size = 16),  # Increase font size for y axis title
    legend.text = element_text(face = "bold", size = 14),  # Increase font size for legend text
    legend.title = element_text(face = "bold", size = 16)  # Increase font size for legend title
  )

# Define the output file path for the Metabolite_Flux plot
output_file_metabolite_flux <- "Metabolite_Flux_Plot.png"

# Save the Metabolite_Flux plot with updated features
ggsave(output_file_metabolite_flux, 
       plot = ggplot(combined_data, aes(x = reorder(HumanReadableName, Value), y = Value, fill = flux_type)) +
         geom_bar(stat = "identity", position = "dodge") +
         coord_flip() +
         labs(
           title = "",
           x = "Metabolites",
           y = "Mean Flux Value",
           fill = "Flux Type"
         ) +
         scale_fill_manual(values = c("Production" = "#4B0082", "Consumption" = "#D3B8E6")) +  # Dark purple and pale purple
         theme_minimal() +
         theme(
           axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 14),  # Increase font size for x axis
           axis.text.y = element_text(face = "bold", size = 14),  # Increase font size for y axis
           axis.title.x = element_text(face = "bold", size = 16),  # Increase font size for x axis title
           axis.title.y = element_text(face = "bold", size = 16),  # Increase font size for y axis title
           legend.text = element_text(face = "bold", size = 14),  # Increase font size for legend text
           legend.title = element_text(face = "bold", size = 16)  # Increase font size for legend title
         ),
       width = 10, height = 7)  # Set desired dimensions for the saved plot

