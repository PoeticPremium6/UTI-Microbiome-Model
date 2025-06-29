library(tidyverse)
library(RColorBrewer)

# Set the path to your data
file_path <- "VF_Expression.csv"

# Read the data
virulence_data <- read_csv(file_path)

# Transform data to long format
virulence_data_long <- virulence_data %>%
  pivot_longer(
    cols = -c(VF, Gene),  # Exclude 'VF' and 'Gene' columns from pivoting
    names_to = "Sample",  # New column to store sample names (e.g., A01, A02)
    values_to = "Expression"  # New column for the expression values
  )

# Filter out rows where Expression is 0
virulence_data_long <- virulence_data_long %>%
  filter(Expression > 0)  # Remove rows where Expression is 0

# Define a custom color palette
custom_palette <- c(
  "#D73027", "#FC8D59", "#FEE090", "#6A4C9C", "#4575B4",
  "#1A9850", "#66BD63", "#A6D96A", "#D9EF8B", "#FFFFBF",
  "#6A4C9C", "#FDAE61", "#F46D43", "#D53E4F", "#9E0142"
)

# Reverse the order of the y-axis by flipping the levels of the 'Gene' factor
virulence_data_long <- virulence_data_long %>%
  mutate(Gene = factor(Gene, levels = rev(levels(factor(Gene)))))

# Create the plot with the reversed y-axis
p <- ggplot(virulence_data_long, aes(x = Sample, y = Gene, color = VF, size = Expression)) +
  geom_point() +
  scale_size_continuous(name = "Expression (FPKM)", range = c(3, 10)) +  # Adjust dot sizes and rename legend
  scale_color_manual(name = "Virulence Factor Category", values = custom_palette) +  # Rename color legend
  theme_minimal() +
  labs(
    # title = "Virulence Factor Expression by Sample",  # Uncomment if a title is needed
    y = "Gene",
    x = "Sample"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),  # Bold and enlarge x-axis ticks
    axis.text.y = element_text(size = 12, face = "bold"),  # Bold and enlarge y-axis ticks
    plot.title = element_text(hjust = 0.5, size = 14),  # Center and style title
    axis.title.x = element_text(size = 14),  # Style x-axis title
    axis.title.y = element_text(size = 14),  # Style y-axis title
    legend.title = element_text(size = 12),  # Style legend title
    legend.text = element_text(size = 10)    # Style legend text
  ) +
  guides(color = guide_legend(override.aes = list(size = 7)))  # Adjust legend dot size for clarity

# Display the plot
print(p)

# Save the plot to the specified directory
output_file_path <- "virulence_factor_expression_ordered_by_VF.png"
ggsave(filename = output_file_path, plot = p, width = 15, height = 10, units = "in")

#################
# Load necessary libraries
library(tidyverse)
library(fgsea)
library(msigdbr)  # For pathway data from MSigDB, including metabolic pathways
library(ggplot2)

# Set the path to your data file
file_path <- "master_gene_abundance_UTI89.csv"

# Read the gene expression matrix
gene_expression <- read_csv(file_path)

# Calculate the average expression per gene across samples as a simple metric
gene_expression <- gene_expression %>%
  rowwise() %>%
  mutate(Average_Expression = mean(c_across(starts_with("A"):starts_with("H")), na.rm = TRUE)) %>%
  ungroup()

# Load necessary libraries
library(dplyr)

# Select the GeneName and Average_Expression columns, then arrange and deframe for fgsea
ranked_genes <- gene_expression %>%
  dplyr::select(GeneName, Average_Expression) %>%
  arrange(desc(Average_Expression)) %>%
  deframe()  # Deframe for fgsea compatibility

# Download metabolic pathway gene sets for E. coli and filter for metabolism
metabolic_gene_sets <- msigdbr(species = "Escherichia coli", category = "C2", subcategory = "CP:KEGG") %>%
  filter(grepl("metabolism", gs_name, ignore.case = TRUE)) %>%
  split(.$gs_name)

# Run fgsea with the updated ranked genes
gsea_results <- fgsea(pathways = metabolic_gene_sets, stats = ranked_genes, nperm = 1000)

# Filter for significant pathways
significant_pathways <- gsea_results %>%
  filter(padj < 0.05) %>%
  arrange(padj)

# Plot top metabolic pathways
plot_data <- significant_pathways %>%
  slice(1:10)  # Select top 10 pathways

# Create a bar plot for top enriched metabolic pathways
p <- ggplot(plot_data, aes(x = reorder(pathway, NES), y = NES)) +
  geom_col(aes(fill = -log10(padj))) +
  coord_flip() +
  scale_fill_gradient(low = "skyblue", high = "darkblue") +
  theme_minimal() +
  labs(title = "Top Metabolic Pathways Enriched in E. coli Gene Expression",
       x = "Metabolic Pathway",
       y = "Normalized Enrichment Score (NES)",
       fill = "-log10 Adjusted p-value") +
  theme(axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# Display the plot
print(p)

# Optionally, save the plot
output_file_path <- "Figure2.png"
ggsave(filename = output_file_path, plot = p, width = 10, height = 8, units = "in")


######################
# Load necessary libraries
library(dplyr)
library(ggplot2)
library(tidyr)

# Load the data
data <- read.delim("master_gene_abundance_UTI89.tsv")

# Clean up the 'Gene_ID' to remove the 'gene-' prefix
data <- data %>%
  mutate(Gene_ID = sub("gene-", "", Gene_ID))

# Create a new column 'Gene_Label' that only replaces '-' in 'Gene_Name' with the cleaned 'Gene_ID'
data <- data %>%
  mutate(Gene_Label = ifelse(Gene_Name == "-", Gene_ID, Gene_Name))

# Convert data to long format
data_long <- data %>%
  pivot_longer(cols = -c(Gene_ID, Gene_Name, Gene_Label), names_to = "Sample", values_to = "Expression")

# Apply log transformation to Expression for better visualization of smaller values
data_long <- data_long %>%
  mutate(Expression = log1p(Expression))  # log1p to handle zeros gracefully

# Group by Sample and arrange genes to find the top-10 unique highly expressed genes for each sample
top_genes <- data_long %>%
  group_by(Sample) %>%
  arrange(desc(Expression)) %>%
  filter(!duplicated(Gene_Label)) %>%  # Keep unique genes based on Gene_Label
  slice_max(Expression, n = 10) %>%  # Select top-10 genes per sample
  ungroup()

# Plotting the top-10 unique highly expressed genes for each sample
# Improved Plot for Top Genes with Enhanced Styling
plot <- ggplot(top_genes, aes(x = reorder(Gene_Label, Expression), y = Expression, fill = Sample)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~ Sample, scales = "free_y") +
  labs(x = "Gene Name or Gene ID", y = "Log-Scaled Expression Level") +
  theme_minimal() +
  theme(
    legend.position = "none",  # Hide legend
    axis.title = element_text(face = "bold", size = 14),  # Bold and larger axis titles
    axis.text.x = element_text(size = 12),  # Larger x-axis tick labels
    axis.text.y = element_text(face = "bold", size = 10),  # Bold and larger y-axis tick labels
    strip.text = element_text(face = "bold", size = 14),  # Bold and larger facet title (subplot headers)
    panel.grid.major.y = element_line(color = "gray", linetype = "dashed"),  # Dashed horizontal gridlines
    panel.grid.major.x = element_blank()  # Remove vertical gridlines
  )


# Save the plot in the specified directory
ggsave("top_transcripts.png", plot, width = 10, height = 8, dpi = 300)

################################
# Load necessary libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(RColorBrewer)

# Set the path to your data
file_path <- "master_subsystem_counts.csv"

# Read the data
heatmap_data <- read.csv(file_path)

# Reshape the data to long format
heatmap_data_long <- heatmap_data %>%
  pivot_longer(cols = -Subsystem, names_to = "Model", values_to = "Count")

# Normalize the counts row-wise by dividing each count by the maximum count in that row
heatmap_data_long <- heatmap_data_long %>%
  group_by(Subsystem) %>%
  mutate(Normalized_Count = Count / max(Count)) %>%
  ungroup()

# Flip the y-axis order by reversing the factor levels for 'Subsystem'
heatmap_data_long$Subsystem <- factor(heatmap_data_long$Subsystem, levels = rev(unique(heatmap_data_long$Subsystem)))

# Set the number of top subsystems to keep
top_n_subsystems <- 15  # Adjust this value as needed

# Filter the data to keep the top N subsystems based on maximum counts
top_subsystems <- heatmap_data_long %>%
  group_by(Subsystem) %>%
  summarize(max_count = max(Count)) %>%
  top_n(top_n_subsystems, max_count) %>%
  pull(Subsystem)

# Filter the original data to keep only the top subsystems
heatmap_data_filtered <- heatmap_data_long %>%
  filter(Subsystem %in% top_subsystems)

# Generate the heatmap with the filtered data
heatmap_plot <- ggplot(heatmap_data_filtered, aes(x = Model, y = Subsystem, fill = Normalized_Count)) +
  geom_tile() +
  scale_fill_gradientn(colors = brewer.pal(9, "Purples"), limits = c(0, 1)) +  # Purple for high values, white for low
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = "bold", vjust = 1),  # Adjust x-axis labels (larger, bold)
    axis.text.y = element_text(size = 12, face = "bold"),  # Adjust y-axis labels (larger, bold)
    axis.title = element_blank(),  # Remove axis titles
    legend.title = element_text(face = "bold"),  # Bold legend title
    legend.text = element_text(size = 10)  # Legend text size
  ) +
  labs(
    fill = "Normalized Count",  # Legend label
    x = "Sample",  # X-axis title
    y = "Metabolic Subsystems"  # Y-axis title
  ) +
  ggtitle("")  # Empty title (can add a title here if desired)

# Display the plot
print(heatmap_plot)

# Save the plot with filtered subsystems in a landscape format
output_file_path <- "heatmap_filtered_subsystems_landscape.png"
ggsave(filename = output_file_path, plot = heatmap_plot, width = 15, height = 8, units = "in")

################################
# Load necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)

# Read in the CSV file
file_path <- "model_summary.csv"
model_data <- read.csv(file_path)

# Check the structure of the data to ensure it has the correct columns
head(model_data)

# Reshape data for plotting
model_data_long <- model_data %>%
  pivot_longer(cols = c(Reactions, Metabolites), names_to = "Type", values_to = "Count")

# Plot using ggplot2
plot <- ggplot(model_data_long, aes(x = ModelIDs, y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") + 
  labs(x = "Model ID", y = "Count", fill = "Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "bold", size = 12),  # Bold x-axis ticks and increased size
        axis.text.y = element_text(face = "bold", size = 12),  # Bold y-axis ticks and increased size
        axis.title.x = element_text(face = "bold", size = 12),  # Bold x-axis title
        axis.title.y = element_text(face = "bold", size = 12),  # Bold y-axis title
        legend.title = element_text(face = "bold"),  # Bold legend title
        legend.text = element_text(face = "bold"),   # Bold legend text
        plot.title = element_blank()) +  # Remove plot title
  scale_fill_manual(values = c("#4B0082", "#9370DB"))  # Dark and light purple

# Save the plot
outputFile <- "/model_comparison_plot.png"
ggsave(outputFile, plot, width = 10, height = 6)

# Print confirmation
print(paste("Plot saved to", outputFile))

######################################
######################################
library(dplyr)
library(tidyr)
library(ggplot2)

# Define the directory path
directory_path <- "/Uromicrobiome/Bacarena"

# List all CSV files in the directory
file_list <- list.files(directory_path, pattern = "*.csv", full.names = TRUE)

# Initialize an empty list to hold the data frames
data_list <- list()

# Loop through each file to read, transpose, and store the data
for (file in file_list) {
  # Read the CSV file
  data <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
  
  # Check if the data frame is not empty and has at least one column
  if (ncol(data) > 0) {
    # Transpose the data so metabolites are in rows and flux values are in columns
    transposed_data <- t(data)  # Transpose the data matrix
    transposed_data <- as.data.frame(transposed_data)  # Convert to a data frame
    
    # Ensure that the transposed data has columns before renaming
    if (ncol(transposed_data) > 0) {
      # Extract metabolites from row names
      metabolites <- rownames(transposed_data)
      
      # Set the column names
      colnames(transposed_data) <- c("Flux_value")  # Name the first column for flux values
      transposed_data$Metabolite <- metabolites  # Add metabolites as a column
      
      # Reorder columns: first Metabolite, then Flux_value
      transposed_data <- transposed_data[, c("Metabolite", "Flux_value")]
      
      # Set the column name based on the file name (without .csv extension)
      sample_name <- gsub(".csv", "", basename(file))
      colnames(transposed_data)[2] <- sample_name  # Rename the flux column to the sample name
      
      # Append to the data list
      data_list[[file]] <- transposed_data
    } else {
      warning(paste("Empty data after transpose in file:", file))
    }
  } else {
    warning(paste("Empty file or malformed data:", file))
  }
}

# Merge all data frames by the "Metabolite" column
if (length(data_list) > 0) {
  merged_data <- Reduce(function(x, y) full_join(x, y, by = "Metabolite"), data_list)
  
  # Optional: Replace any NA values with zeros or another placeholder if required
  merged_data[is.na(merged_data)] <- 0
  
  # Print the merged data (first few rows)
  print(head(merged_data))
  
  # Optionally, save the merged data to a CSV file
  write.csv(merged_data, "merged_flux_data.csv", row.names = FALSE)
} else {
  print("No valid data found to merge.")
}

# Control column for comparison
control_column <- "Control_Flux_UTI89"

# Function to calculate the change from control, perform statistical test, and extract the top 5 positive and negative metabolites for each sample
compare_metabolite_change_per_sample <- function(data, control_column) {
  # Extract metabolite names (rows)
  metabolites <- data$Metabolite
  
  # Subset the data (excluding the 'Metabolite' column)
  data_sub <- data %>% select(-Metabolite)
  
  # Ensure the control column is excluded from the difference calculations for itself
  data_sub <- data_sub[, colnames(data_sub) != control_column]
  
  # Create an empty data frame to store results
  top_results <- data.frame(Sample = character(),
                            Metabolite = character(),
                            Change = numeric(),
                            Type = character(),
                            P_value = numeric())
  
  # Iterate through each sample column
  for (sample_name in colnames(data_sub)) {
    # Extract the metabolite values for the current sample and the control
    sample_values <- data[[sample_name]]
    control_values <- data[[control_column]]
    
    # Ensure that both sample_values and control_values are numeric and handle missing data
    if (all(is.numeric(sample_values), is.numeric(control_values))) {
      # Perform t-test for each metabolite (compare sample and control as entire vectors)
      p_values <- sapply(1:nrow(data), function(i) {
        test_result <- tryCatch(
          {
            t.test(sample_values, control_values)$p.value
          },
          error = function(e) {
            return(NA)  # Return NA if t-test fails
          }
        )
        return(test_result)
      })
      
      # Calculate the difference between the sample and the control (can use difference or ratio)
      difference_data <- data_sub[[sample_name]] - data[[control_column]]
      
      # Combine the metabolites, their changes, and p-values into a data frame
      diff_df <- data.frame(Metabolite = metabolites, Change = difference_data, P_value = p_values)
      
      # Sort by absolute change
      sorted_diff <- diff_df %>% arrange(desc(abs(Change)))
      
      # Get top 5 positive and negative changes
      top_5_positive <- sorted_diff %>% filter(Change > 0) %>% head(2)
      top_5_negative <- sorted_diff %>% filter(Change < 0) %>% head(2)
      
      # Add 'Sample' and 'Type' columns to indicate positive or negative change
      top_5_positive$Sample <- sample_name
      top_5_positive$Type <- "Positive"
      top_5_negative$Sample <- sample_name
      top_5_negative$Type <- "Negative"
      
      # Combine positive and negative results
      sample_results <- rbind(top_5_positive, top_5_negative)
      
      # Append to the final results dataframe
      top_results <- rbind(top_results, sample_results)
    } else {
      warning(paste("Non-numeric values found in sample or control for sample:", sample_name))
    }
  }
  
  # Return the final table with top 5 positive and negative metabolites for each sample and their p-values
  return(top_results)
}

# Apply the function to compare each sample to the control and include statistical test results
comparison_results_per_sample <- compare_metabolite_change_per_sample(merged_data, control_column)

# Print the comparison results
print(comparison_results_per_sample)

# Optionally, write the results to a CSV
write.csv(comparison_results_per_sample, "/Bacarena/top_metabolites_with_p_values_per_sample.csv", row.names = FALSE)

##########
library(dplyr)
library(ggplot2)

# Load the mapping file to translate metabolite IDs to human-readable names
mapping_file <- read.csv("metabolites_mapping.csv")

# Clean the Human_Readable_Name column in the mapping file (remove any special characters or extra spaces)
mapping_file$Human_Readable_Name <- trimws(mapping_file$Human_Readable_Name)

# If 'Human_Readable_Name' was renamed, we will fix it here
if ("Human_Readable_Name.x" %in% colnames(top_fluxes_per_sample)) {
  top_fluxes_per_sample <- top_fluxes_per_sample %>%
    dplyr::rename(Human_Readable_Name = "Human_Readable_Name.x")  # Rename the column back to 'Human_Readable_Name'
}

# Merge comparison results with the mapping file
top_fluxes_per_sample <- comparison_results_per_sample %>%
  left_join(mapping_file, by = c("Metabolite" = "Metabolite_ID"))

# Print the column names for confirmation
print(colnames(top_fluxes_per_sample))

# Remove any column that starts with "Human_Readable_Name." or other unwanted columns (if necessary)
top_fluxes_per_sample <- top_fluxes_per_sample %>%
  dplyr::select(-starts_with("Human_Readable_Name."))  # Remove any 'Human_Readable_Name.*' columns

# Check if 'Change' column exists, then calculate log10_Change
if ("Change" %in% colnames(top_fluxes_per_sample)) {
  # Add log10 transformation for change values
  top_fluxes_per_sample <- top_fluxes_per_sample %>%
    mutate(log10_Change = ifelse(Change > 0, log10(Change), -log10(abs(Change))))  # Log-transform Change
} else {
  warning("The 'Change' column is missing from the dataset.")
}

# Assign 'Type' based on whether the 'Change' is positive or negative
top_fluxes_per_sample$Type <- ifelse(top_fluxes_per_sample$Change > 0, "Produce", "Consume")

# Print the cleaned data with the new log10_Change and Type columns
head(top_fluxes_per_sample)

# Plot the top 2 positive and negative changes per sample with log10 transformation
plot <- ggplot(top_fluxes_per_sample, aes(x = reorder(Human_Readable_Name, log10_Change), y = log10_Change, fill = Type)) +
  geom_bar(stat = "identity", show.legend = TRUE) +
  geom_text(aes(label = sprintf("p = %.2e", P_value)), 
            position = position_stack(vjust = 0.5),  # Center the p-value labels within the bars
            size = 2.5,  # Reduce the font size for p-values
            color = "white") +  # Make p-values white for better contrast
  coord_flip() +  # Flip for better readability of metabolite names
  facet_wrap(~ Sample, scales = "free_y") +  # Facet by simplified sample
  labs(
    title = "Top 2 Metabolites Produced and Consumed per Sample (Log10 Scale)",
    x = "Metabolite",
    y = "Log10 Change (Flux)"
  ) +
  scale_fill_manual(values = c("Produce" = "forestgreen", "Consume" = "darkorange")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  # Rotate and bold x-axis ticks
    axis.text.y = element_text(face = "bold"),  # Bold y-axis ticks
    axis.title.x = element_text(face = "bold", size = 12),  # Bold x-axis label
    axis.title.y = element_text(face = "bold", size = 10),  # Bold y-axis label
    #plot.title = element_text(face = "bold", size = 14),  # Bold title
    strip.text = element_text(face = "bold", size = 10)  # Bold facet labels
  )

# Display the plot
print(plot)

# Optionally, save the plot
ggsave("top_fluxes_per_sample_plot.png", plot, width = 10, height = 8)
