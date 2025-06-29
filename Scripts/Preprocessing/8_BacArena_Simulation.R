# Import required libraries
library("BacArena")
library("R.matlab")
library("parallel")

# Set working directory to where your models and diet information are stored
setwd("/rUTI/Context-Specific_Models")

# Function to read in a MATLAB model and prepare it for BacArena
prepareModel <- function(model_name) {
  model_path <- paste0(model_name, "_context.mat")
  model <- readMATmod(model_path)
  Bac(model, setAllExInf = TRUE, cellweight_sd = 0)
}

# Load diet information
diet <- read.table("/rUTI/urine_media.csv", sep = ",", header = TRUE)

# Prepare models
model_names <- c("Clostridium_cocleatum", "Collinsella_intestinalis", "Cronobacter_turicensis",
                 "Desulfomicrobium_orale", "Escherichia_coli", "Kosakonia_sacchari",
                 "Neobacillus_niacini", "Odoribacter_laneus", "Pantoea_agglomerans",
                 "Salmonella_enterica", "Serratia_marcescens", "Streptococcus_troglodytae",
                 "Yersinia_enterocolitica", "Truepera_radiovictrix", "Prevotella_oris",
                 "Slackia_faecicanis", "Slackia_piriformis")

models <- lapply(model_names, prepareModel)

# Set up the simulation environment (Arena)
arena <- Arena(n = 100, m = 100)
for (i in 1:length(models)) {
  arena <- addOrg(arena, models[[i]], amount = 10) # Adjust amount as per your specific study
}

# Add substances from diet to the Arena
arena <- addSubs(arena, mediac = paste0("EX_", diet$compounds, "_e0"), smax = diet$maxFlux, unit = "mM")

# Run the simulation
sim <- simEnv(arena, time = 4) # Adjust time as needed

# Save the simulation results
saveRDS(sim, file = "/rUTI/Community_Model/Simulation_Results.rds")

# Analysis and plotting
plotAbundance(sim) + ggplot2::scale_fill_manual(values = colpal2) + ggplot2::scale_color_manual(values = colpal2)
plotGrowthCurve(sim)[[1]]
plotSubCurve(sim)[[1]]
findFeeding3(sim, time = 4, mets = names(getVarSubs(sim)))

# Export variable substances data to CSV
sim <- readRDS(file = "/rUTI/Community_Model/Simulation_Results.rds")
df <- getVarSubs(sim)
write.csv(df, file = '/rUTI/Community_Model/Simulation_Results.csv')
