#!/bin/bash

# OT Mapper: Map Off-Target (OT) candidates to gene annotations
# Usage: ./script_name.sh input_bed_file

# Check if input BED file is provided
BED1=${1?Error: No input BED file provided. Usage: ./script_name.sh input_bed_file}

# Print start message
echo "Mapping $BED1 with GENCODE V39 exon and intron annotations..."

# Create timestamp for the analysis
timestamp=$(date +"%Y%m%d_%H%M%S")
echo "Analysis started at: $timestamp"

# Create output directory
mkdir -p analysis/OT_mapper_results/

# Perform intersection with exon and intron databases
# -a: Input BED file (OT candidates)
# -b: Reference BED files (exons and introns)
# -wb: Write the original entry in B for each overlap
# The output is then processed to keep only relevant fields
echo "Intersecting OT candidates with exon and intron annotations..."
bedtools intersect \
    -a $BED1 \
    -b data/UCSC_exons_modif_canonical.bed data/UCSC_introns_modif_canonical.bed \
    -wb | cut -f 5,1-3,8-9 > analysis/OT_mapper_results/OT_mapped.tsv

# Print completion message
echo "Mapping completed."
echo "Results saved in: analysis/OT_mapper_results/OT_mapped.tsv"

# Create end timestamp
end_timestamp=$(date +"%Y%m%d_%H%M%S")
echo "Analysis finished at: $end_timestamp"

echo "Done!"