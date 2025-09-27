# Load the necessary libraries
library(optparse)

# Define the command line options
option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Input file name"),
  make_option(c("-o", "--output"), type="character", help="Output file name")
)

# Parse the command line options
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Read the input file into a data frame
table <- read.table(opt$input, header=FALSE, sep="\t")

# Load the function
convert_chr_to_refseq <- function(chr, ref_table) {
  refseq_id <- ref_table[match(chr, ref_table[, 1]), 2]
  if (is.na(refseq_id)) {
    stop(paste0("Invalid chromosome identifier: ", chr))
  }
  return(refseq_id)
}

ref_table <- read.table(text="
chr1	NC_000001.11
chr2	NC_000002.12
chr3	NC_000003.12
chr4	NC_000004.12
chr5	NC_000005.10
chr6	NC_000006.12
chr7	NC_000007.14
chr8	NC_000008.11
chr9	NC_000009.12
chr10	NC_000010.11
chr11	NC_000011.10
chr12	NC_000012.12
chr13	NC_000013.11
chr14	NC_000014.9
chr15	NC_000015.10
chr16	NC_000016.10
chr17	NC_000017.11
chr18	NC_000018.10
chr19	NC_000019.10
chr20	NC_000020.11
chr21	NC_000021.9
chr22	NC_000022.11
chrX	NC_000023.11
chrY	NC_000024.10", header = FALSE, stringsAsFactors = FALSE)

# Apply the function to the first column of the data frame to get RefSeq IDs
table$V1 <- sapply(table$V1, convert_chr_to_refseq, ref_table = ref_table)

# Construct the Primer-BLAST URL for each row of the table
urls <- character(nrow(table))
for (i in 1:nrow(table)) {
    url <- paste0("https://www.ncbi.nlm.nih.gov/tools/primer-blast/index.cgi?LINK_LOC=bookmark&INPUT_SEQUENCE=", table$V1[i], "&PRIMER5_START=", table$V2[i]-500, "&PRIMER5_END=", table$V2[i]-150, "&PRIMER3_START=", table$V3[i]+150, "&PRIMER3_END=", table$V3[i]+500, "&OVERLAP_5END=7&OVERLAP_3END=4&PRIMER_PRODUCT_MIN=400&PRIMER_PRODUCT_MAX=800&PRIMER_NUM_RETURN=10&PRIMER_MIN_TM=56&PRIMER_OPT_TM=58&PRIMER_MAX_TM=62&PRIMER_MAX_DIFF_TM=2&PRIMER_ON_SPLICE_SITE=0&SEARCHMODE=0&SPLICE_SITE_OVERLAP_5END=7&SPLICE_SITE_OVERLAP_3END=4&SPLICE_SITE_OVERLAP_3END_MAX=8&SPAN_INTRON=off&MIN_INTRON_SIZE=1000&MAX_INTRON_SIZE=1000000&SEARCH_SPECIFIC_PRIMER=on&EXCLUDE_ENV=on&EXCLUDE_XM=off&TH_OLOGO_ALIGNMENT=off&TH_TEMPLATE_ALIGNMENT=off&ORGANISM=Homo%20sapiens&PRIMER_SPECIFICITY_DATABASE=PRIMERDB/genome_selected_species&TOTAL_PRIMER_SPECIFICITY_MISMATCH=1&PRIMER_3END_SPECIFICITY_MISMATCH=1&MISMATCH_REGION_LENGTH=5&TOTAL_MISMATCH_IGNORE=6&MAX_TARGET_SIZE=4000&ALLOW_TRANSCRIPT_VARIANTS=off&HITSIZE=50000&EVALUE=30000&WORD_SIZE=7&MAX_CANDIDATE_PRIMER=500&PRIMER_MIN_SIZE=18&PRIMER_OPT_SIZE=22&PRIMER_MAX_SIZE=25&PRIMER_MIN_GC=20.0&PRIMER_MAX_GC=80.0&GC_CLAMP=1&NUM_TARGETS_WITH_PRIMERS=1000&NUM_TARGETS=20&MAX_TARGET_PER_TEMPLATE=100&POLYX=5&SELF_ANY=8.00&SELF_END=3.00&PRIMER_MAX_END_STABILITY=9&PRIMER_MAX_END_GC=2&PRIMER_MAX_TEMPLATE_MISPRIMING_TH=40.00&PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH=70.00&PRIMER_MAX_SELF_ANY_TH=45.0&PRIMER_MAX_SELF_END_TH=35.0&PRIMER_PAIR_MAX_COMPL_ANY_TH=45.0&PRIMER_PAIR_MAX_COMPL_END_TH=35.0&PRIMER_MAX_HAIRPIN_TH=24.0&PRIMER_MAX_TEMPLATE_MISPRIMING=12.00&PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING=24.00&PRIMER_PAIR_MAX_COMPL_ANY=8.00&PRIMER_PAIR_MAX_COMPL_END=3.00&PRIMER_MISPRIMING_LIBRARY=repeat/repeat_9606&NO_SNP=on&LOW_COMPLEXITY_FILTER=on&MONO_CATIONS=50.0&DIVA_CATIONS=1.5&CON_ANEAL_OLIGO=50.0&CON_DNTPS=0.6&SALT_FORMULAR=1&TM_METHOD=1&PRIMER_INTERNAL_OLIGO_MIN_SIZE=18&PRIMER_INTERNAL_OLIGO_OPT_SIZE=20&PRIMER_INTERNAL_OLIGO_MAX_SIZE=27&PRIMER_INTERNAL_OLIGO_MIN_TM=57.0&PRIMER_INTERNAL_OLIGO_OPT_TM=60.0&PRIMER_INTERNAL_OLIGO_MAX_TM=63.0&PRIMER_INTERNAL_OLIGO_MAX_GC=80.0&PRIMER_INTERNAL_OLIGO_OPT_GC_PERCENT=50&PRIMER_INTERNAL_OLIGO_MIN_GC=20.0&PICK_HYB_PROBE=off&NEWWIN=on&NEWWIN=on&SHOW_SVIEWER=true")
    urls[i] <- url
}

# Add the URLs to the output table
table$V6 <- urls

# Write the output table to a file
write.table(table, opt$output, row.names=FALSE, col.names=FALSE, sep="\t")
