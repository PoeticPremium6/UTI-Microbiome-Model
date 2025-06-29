# Metatranscriptomics-based metabolic modeling of patient-specific urinary microbiome during infection
Leveraging patient-specific metatranscriptomes to reconstruct and simulate microbiome community model

Our preprint is available here:
https://www.biorxiv.org/content/10.1101/2024.03.25.586446v1.full.pdf

![Screenshot 2024-09-18 225507](https://github.com/user-attachments/assets/0d2284c8-6f6d-4d0d-bc2b-8e374c1011db)



# UTI Microbiome Metabolic Modeling Pipeline
# Intro
This repository implements a metatranscriptomics-to-metabolic-modeling pipeline for patient urinary microbiomes, as described in the manuscript. Raw RNA sequencing reads are first processed with standard bioinformatics tools: FastQC for read quality assessment, Cutadapt for adapter trimming, and PRINSEQ-lite for filtering and quality trimming. Host (human) reads are removed by aligning against the GRCh38 reference with Bowtie2. Ribosomal RNA and mRNA reads are then identified and filtered using SortMeRNA with the SILVA and Rfam databases. This yields high-quality, non-human mRNA and rRNA reads for downstream analysis. 

# How the Repository Runs
This repository is organized into three main script directories:
1) Preprocessing: Contains scripts for raw data processing, including quality control, rRNA/mRNA separation , mapping and expression quantification, and metabolic model reconstruction. These scripts prepare the inputs for downstream metabolic modeling and community simulations.
2) Main_Analysis: Contains R scripts (Figure2.R–Figure6.R) used to generate the main manuscript figures from the processed data and simulation outputs. Each script is self-contained and annotated with input/output file paths and statistical methods.
3) Supplementary_Analysis: Includes scripts and data used to reproduce supplementary figures or perform additional exploratory analyses that support the main text.


# rRNA and mRNA Analysis
Taxonomic profiling of both metatranscriptomic rRNA and 16S rDNA data is performed by clustering sequences into OTUs using USEARCH and classifying them with against the RDP 16S database. The RDP classifier is a Naïve Bayesian method for rapid 16S taxonomy assignment. The R package phylose is used to merge OTU tables with sample metadata, filter low-abundance taxa, and normalize counts. Gene expression quantification proceeds by mapping the filtered mRNA reads to species reference genomes (see Supplementary Table S3). We build HISAT2 genome indexes and align reads with HISAT2 in RNA-seq mode. The alignments are processed with SAMtools. Gene-level counts are obtained with FeatureCounts. StringTie is then used to assemble transcripts and calculate expression values (FPKM) for each gene. All species’ gene counts are combined into cohort-wide tables. In particular, E. coli counts are compared against the UTI89 genome to identify expressed virulence genes using the Virulence Factor Database (VFDB). 

# Context-specific Metabolic Modelling
For metabolic modeling, we build species-specific genome-scale models using gapseq with each bacterial genome from NCBI. We perform gapseq’s “find” and “find-transport” steps to predict reactions and transporters, then create draft models and gap-fill against our custom urine medium. The custom “in silico urine” medium is based on the Human Urine Metabolome (Bouatra et al. 2013) which identified 445 urinary metabolites. Metabolite concentrations were normalized to creatinine levels and formatted as VMH exchange reactions. We then integrate expression data into the GEMs using the COBRA Toolbox  in MATLAB. Gene expression (FPKM) is mapped to model reactions, and the fast-consistency algorithm (fastcc) is applied to identify a flux-consistent subnetwork. Next, the fastcore algorithm extracts a minimal “active” reaction set: reactions with expression in the top 25% are treated as core, and fastcore expands this to a flux-consistent model. The resulting context-specific models are  tested by flux balance analysis and flux variabilility analysis. Finally, we simulate patient-specific microbial communities with BacArena in R. Community models are assembled by combining each patient’s species-specific GEMs in a 100×100  arena with initial abundances proportional to the OTU counts. The urine medium (VMH-based) provides extracellular metabolites. Simulations run for several time steps (1-hour increments, 4 hours total, in triplicate) to capture growth and metabolic interactions. We analyze outputs for species growth, metabolite production/consumption, and cross-feeding patterns between species.

# Software Dependencies
The pipeline requires the following software:
FastQC v0.11.8
Cutadapt v1.5
PRINSEQ-lite v0.20.4
Bowtie2 v2.5.4
SortMeRNA v4.3.7 
CD-HIT v4.8.1
USEARCH v12.0
HISAT2 v2.1.0 
SAMtools v1.16
featureCounts (Subread) v2.0.3 
StringTie v2.1.4
Python v3.x
gapseq v1.3.1
MATLAB (R2019 with COBRA Toolbox (v3.0)
R v4.3.2  
Phyloseq (v1.5)
BacArena (v1.8.1)
(Other R packages such as dplyr, ggplot2, etc., are also used in analysis scripts.)
The pipeline was developed on Linux; Unix-compatible environments (macOS/Linux) are recommended. Some tools (Cutadapt, gapseq, etc.) have their own dependencies (e.g. Python libraries, BLAST for gapseq). Users should ensure the above software is installed and in $PATH.

# Reference Databases
The following reference data are required:
Human genome (GRCh38) 
Rfam database
RDP 16S database 
NCBI Genomes and Annotations
Virulence Factor Database (VFDB)
Human Urine Metabolome (HMDB)
Virtual Metabolic Human (VMH)
