#!/bin/bash

# Define directories
input_dir="/mRNA"
output_dir="/CD-HIT"
cdhit_path="/Software/cd-hit-v4.8.1/cd-hit"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Define sample names
sample_names=("A01" "A02" "B01" "B02" "C01" "C02" "D01" "D02"
              "E01" "E02" "F01" "F02" "G01" "H01" "H25361"
              "H25362" "H25363" "H25364" "H25365")

# Loop through each sample and run CD-HIT with increased memory
for sample in "${sample_names[@]}"; do
    input_file="${input_dir}/${sample}_mRNA.fastq"
    output_file="${output_dir}/${sample}_mRNA_cdhit"
    log_file="${output_dir}/${sample}_mRNA_cdhit.log"

    # Run CD-HIT and log the output with increased memory
    "$cdhit_path" -i "$input_file" -o "$output_file" -c 0.95 -n 5 -M 4000 > "$log_file" 2>&1

    echo "CD-HIT completed for $sample. Results saved in $output_file, log saved in $log_file."
done

