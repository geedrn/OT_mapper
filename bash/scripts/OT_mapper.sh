#!/bin/bash

# Load OT candidates
BED1=${1?Error: no input BED}

echo Mapping $BED1 with GENCODE V39 exon and intron annotations...
echo ...
# Create time stamp
echo ...
# Intersect with two databases
mkdir -p analysis/OT_mapper_results/
bedtools intersect -a $BED1 -b data/UCSC_exons_modif_canonical.bed data/UCSC_introns_modif_canonical.bed  -wb | cut -f 5,1-3,8-9 > analysis/OT_mapper_results/OT_mapped.tsv
echo ...
echo Done!
