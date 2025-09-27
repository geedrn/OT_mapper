#!/bin/sh

# OfferOT v1.0 - Off-target analysis script for CRISPR guide RNAs
# Created by RN

# Function to download and process data for a given strand
process_strand() {
    local strand=$1
    local sp=$2
    local seed=$3
    local suffix=$4

    # Download data from GGGenome
    wget -q "https://gggenome.dbcls.jp/ja/hg38/3/${strand}/nogap/${sp}.csv" -O "${suffix}_1.csv"
    wget -q "https://gggenome.dbcls.jp/ja/hg38/1/${strand}/nogap/${seed}.csv" -O "${suffix}_2.csv"

    echo "${suffix}_1.csv and ${suffix}_2.csv were imported from GGGenome."

    # Process downloaded files
    for i in 1 2; do
        cut -f 9 -d "," "${suffix}_${i}.csv" > "OT${i}_${suffix}.csv"
        sed '1,5d' "OT${i}_${suffix}.csv" | sed 's/"//g' > "OT${i}_processed_${suffix}.csv"
    done

    # Find alignments
    grep -f "OT2_processed_${suffix}.csv" "OT1_processed_${suffix}.csv" > "OT_alignment_${suffix}.csv"

    # Process alignments based on strand
    if [ "$strand" = "+" ]; then
        grep "AG$" "OT_alignment_${suffix}.csv" > "AG.csv"
        grep "GG$" "OT_alignment_${suffix}.csv" > "GG.csv"
        cat AG.csv GG.csv > "${suffix}.csv"
    else
        grep "^CC" "OT_alignment_${suffix}.csv" > "CC.csv"
        grep "^CT" "OT_alignment_${suffix}.csv" > "CT.csv"
        cat CC.csv CT.csv > "${suffix}.csv"
    fi
}

# Prompt user for gRNA sequence (spacer)
read -p 'Input gRNA sequence (DO NOT INCLUDE PAM): ' Spacer

# Validate spacer length
if [ ${#Spacer} != 20 ]; then
    echo "Error: Input needs to be 20 nt!"
    exit 1
fi

echo "This is OfferOT v1.0 created by RN"
echo "Spacer: $Spacer is imported successfully"

# Prompt user for seed sequence length
read -p 'Choose the length of seed sequence (8 or 12 nt): ' VAR

# Validate seed sequence length
if [ "$VAR" != "8" ] && [ "$VAR" != "12" ]; then
    echo "Error: Input needs to be 8 or 12!"
    exit 1
fi

# Prepare spacer and PAM sequence
SandP="${Spacer}NGG"

# Extract seed sequence based on user input
if [ "$VAR" = "12" ]; then
    Seed=$(echo $SandP | cut -c 9-23)
else
    Seed=$(echo $SandP | cut -c 13-23)
fi

# Process data for both strands
process_strand "+" "$SandP" "$Seed" "plus"
process_strand "-" "$SandP" "$Seed" "minus"

echo "The candidates of OT were imported successfully."

# Display detected off-targets
echo "OTs were detected successfully."
echo "=============="
cat plus.csv minus.csv
echo "=============="

# Combine and process final results
cat plus.csv minus.csv > combined.csv
grep -f combined.csv plus_1.csv minus_1.csv > OT_list_final.csv

# Generate UCSC compatible list and BED file
awk -F ',' '{print $9 "," $1 ":" $6 "-" $7}' OT_list_final.csv | sed -e 's/"//g' > UCSC_list_final.csv
awk -F ',' '{print $1"\t"$3"\t"$4}' OT_list_final.csv | sed -e 's/"//g' > OT_candidate.bed

# Create directories and move files
mkdir -p analysis/OT analysis/intermediate
mv OT_list_final.csv UCSC_list_final.csv analysis/OT/
mv *.csv analysis/intermediate/
mv OT_candidate.bed analysis/OT/