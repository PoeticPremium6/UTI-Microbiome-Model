conda activate gapseq-env

# Find transporters and reactions
gapseq find -p all -b 100 -t 8 /rUTI/Genomes/UTI89_genome.fasta

# Reconstruct the metabolic model
gapseq draft -r /rUTI/Gapseq/gapseq_result/UTI89_reactions.tbl -t gapseq_result/UTI89_transporters.tbl -c /rUTI/urine_media -o gapseq_model/UTI89
