# Microbial Community Analysis Script with ENS Calculation

# Load Required Libraries
library(phyloseq)
library(ggplot2)
library(dplyr)
library(vegan)
library(ggrepel)
library(patchwork)

# Define Paths to Data
otu_path <- "/rUTI/otu_table.tsv"
tax_path <- "/rUTI/taxonomy_table.tsv"

# Load and Preprocess Data
otu <- read.delim(otu_path, check.names = FALSE, sep = "\t")
tax <- read.delim(tax_path, check.names = FALSE, sep = "\t", row.names = 1)

# Prepare Phyloseq Object
sample_data <- data.frame(SampleID = rownames(otu), Row.names = rownames(otu))
physeq <- phyloseq(otu_table(as.matrix(otu), taxa_are_rows = TRUE),
                   tax_table(as.matrix(tax)),
                   sample_data(sample_data))

# Beta Diversity and ENS Calculation
physeq_ord <- ordinate(physeq, method = "PCoA", distance = "bray")

# Shannon Diversity and ENS
shannon_div <- estimate_richness(physeq, measures = "Shannon")$Shannon
ens <- exp(shannon_div)  # Exponential of Shannon's index for ENS

# Visualization and Analysis
# Beta Diversity Plot
plot_ordination(physeq, physeq_ord, color = "SampleID") + geom_point()

# Effective Number of Species (ENS) Distribution
ens_plot <- ggplot(data.frame(ENS = ens, SampleID = names(ens)), aes(x = SampleID, y = ENS)) +
            geom_bar(stat = "identity") + theme_minimal() + coord_flip() +
            labs(x = "Sample", y = "Effective Number of Species (ENS)")
print(ens_plot)

# Enterobacteriaceae Diversity vs. PCA Axis 1 (including ENS context)
physeq_entero <- subset_taxa(physeq, Family == "Enterobacteriaceae")
entero_ens <- exp(estimate_richness(physeq_entero, measures = "Shannon")$Shannon)
pca_scores <- scores(physeq_ord)$sites
pca_vs_entero_ens <- ggplot(aes(x = pca_scores[,1], y = entero_ens)) +
                     geom_point() + geom_smooth(method = "lm") +
                     labs(x = "PCA Axis 1", y = "Enterobacteriaceae ENS")
print(pca_vs_entero_ens)

# Save Plots and Results
ggsave("/rUTI/ens_plot.png", ens_plot, width = 10, height = 8, dpi = 300)
ggsave("/rUTI/pca_vs_entero_ens.png", pca_vs_entero_ens, width = 10, height = 8, dpi = 300)
