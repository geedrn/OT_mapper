# OT Mapper

A tool for identifying CRISPR-Cas9 off-target (OT) candidates and annotating them with gene information (exon/intron). This software provides both **Bash** and **R** implementations for flexible usage.

## Overview

**Purpose**: Identify off-target candidates for CRISPR-Cas9 gRNA (SpCas9/NGG), annotate them with gene information (exon/intron), and generate Primer-BLAST links for primer design.

**Key Features**:
- Detects off-target candidates using GGGenome (hg38)
- Filters candidates based on seed sequence overlap and PAM matching
- Annotates candidates with exon/intron information from GENCODE V39
- Generates Primer-BLAST URLs for experimental validation

## Software Concept and Algorithm

### Core Algorithm

The OT Mapper uses a **two-stage filtering approach** to identify high-confidence off-target sites:

1. **Broad Search with Seed Overlap Filtering**
   - Searches the human genome (hg38) for sequences matching the full gRNA spacer+PAM (23nt) allowing up to **3 mismatches**
   - Simultaneously searches for the seed sequence+PAM (8-12nt + PAM) allowing up to **1 mismatch**
   - Uses **coordinate-based overlap detection** (via `bedtools intersect`) to identify full-length candidates that contain a matching seed sequence at the same genomic location
   - This ensures that off-target sites have high similarity in the critical seed region (3' end of the spacer)

2. **PAM Exact Match Filtering**
   - Filters candidates to retain only those with **exact PAM matches** (NGG for SpCas9)
   - For plus strand: sequences ending with AG or GG
   - For minus strand: sequences starting with CC or CT (reverse complement of NGG)

### Why This Approach?

**Biological Rationale**:
- The **seed region** (3' end of the gRNA, 8-12 nucleotides) is critical for Cas9 binding specificity
- Off-target sites with mismatches in the seed region are less likely to be cleaved
- However, sites with seed matches but mismatches elsewhere can still be cleaved
- The **PAM sequence** (NGG) is absolutely required for Cas9 recognition

**Computational Efficiency**:
- Searching for full 23nt sequences with 3 mismatches yields many candidates
- Searching for shorter seed sequences (8-12nt) with 1 mismatch is more specific
- Overlap detection identifies candidates that satisfy both criteria simultaneously
- This reduces false positives while maintaining sensitivity

### Workflow Details

```
Input: 20nt gRNA spacer (e.g., "GCTGAAGCACTGCACGCCGT")
       ↓
Append PAM: "GCTGAAGCACTGCACGCCGTNGG" (23nt)
Extract seed: Last 8-12nt + PAM (e.g., "GCACGCCGTNGG" for 12nt seed)
       ↓
┌─────────────────────────────────────────────────────────┐
│ Stage 1: GGGenome Search                                │
├─────────────────────────────────────────────────────────┤
│ Full sequence (23nt): Search with ≤3 mismatches         │
│ Seed sequence (11-15nt): Search with ≤1 mismatch        │
│ Both searches performed on + and - strands separately    │
└─────────────────────────────────────────────────────────┘
       ↓
┌─────────────────────────────────────────────────────────┐
│ Stage 2: Coordinate-Based Overlap Detection             │
├─────────────────────────────────────────────────────────┤
│ Convert results to BED format (chrom, start, end)       │
│ Use bedtools intersect to find:                          │
│   - Full sequence candidates that overlap with           │
│     seed sequence candidates at same genomic location   │
│ This ensures seed region matches within full candidates │
└─────────────────────────────────────────────────────────┘
       ↓
┌─────────────────────────────────────────────────────────┐
│ Stage 3: PAM Exact Match Filtering                      │
├─────────────────────────────────────────────────────────┤
│ Filter sequences to keep only exact PAM matches:        │
│   Plus strand: End with AG or GG                        │
│   Minus strand: Start with CC or CT                     │
└─────────────────────────────────────────────────────────┘
       ↓
┌─────────────────────────────────────────────────────────┐
│ Stage 4: Gene Annotation                                │
├─────────────────────────────────────────────────────────┤
│ Map candidates to GENCODE V39 annotations using         │
│ bedtools intersect with exon/intron BED files           │
│ Identify which candidates are in exons vs introns       │
└─────────────────────────────────────────────────────────┘
       ↓
┌─────────────────────────────────────────────────────────┐
│ Stage 5: Primer Design Links                            │
├─────────────────────────────────────────────────────────┤
│ Generate Primer-BLAST URLs for each candidate           │
│ Convert chromosome coordinates to RefSeq IDs            │
│ Create links for experimental validation                │
└─────────────────────────────────────────────────────────┘
       ↓
Output: Annotated off-target candidates with primer links
```

### Key Differences from Simple Sequence Matching

**Traditional approach** (string matching):
- Searches for sequences matching the gRNA
- May miss candidates with seed matches but other mismatches
- May include false positives with sequence similarity but wrong genomic context

**OT Mapper approach** (coordinate-based overlap):
- Ensures seed region matches are within full-length candidates
- Uses genomic coordinates to verify spatial relationships
- More accurate and biologically relevant results
- Both Bash and R versions now use this approach for consistency

## Definitions and Parameters

### Core Definitions

**Spacer (gRNA spacer)**:
- The 20-nucleotide guide RNA sequence that targets a specific genomic location
- Does NOT include the PAM sequence
- Standard length: **20 nucleotides** (for SpCas9)
- Example: `GCTGAAGCACTGCACGCCGT`

**PAM (Protospacer Adjacent Motif)**:
- A short DNA sequence required for Cas9 recognition
- Default: **NGG** (any nucleotide followed by two guanines) for SpCas9
- Position: Immediately 3' of the spacer sequence
- The tool automatically appends PAM to the spacer for searching

**Seed Sequence**:
- The 3' end region of the spacer that is critical for Cas9 binding specificity
- Options: **8 to 12 nucleotides** (user-selectable, inclusive)
- Position: Last 8-12 nt of the spacer + PAM
- For 12nt seed: positions 9-23 of spacer+PAM (last 12nt of spacer + 3nt PAM)
- For 8nt seed: positions 13-23 of spacer+PAM (last 8nt of spacer + 3nt PAM)
- For intermediate lengths (9-11nt): positions calculated accordingly
- Mismatches in the seed region significantly reduce off-target activity

**Full Sequence**:
- Complete spacer + PAM sequence (typically 23 nucleotides for 20nt spacer + 3nt PAM)
- Used for broad search with higher mismatch tolerance

**Off-Target (OT) Candidate**:
- A genomic location that matches the gRNA sequence with some mismatches
- Must satisfy both seed overlap and PAM matching criteria

### Search Parameters

**Genome Assembly**:
- Default: **hg38** (human genome reference GRCh38)
- GGGenome API supports: hg38, hg19, mm10, etc.
- Can be changed via `-g` or `--genome` option

**Mismatch Tolerance**:
- **Full sequence mismatch**: Default **3 mismatches**
  - Allows up to 3 nucleotide differences in the full spacer+PAM sequence
  - Searches for sequences with high similarity but not perfect matches
  - GGGenome API parameter: `/3/` in URL
  
- **Seed sequence mismatch**: Default **1 mismatch**
  - Allows up to 1 nucleotide difference in the seed+PAM sequence
  - More stringent requirement for the critical seed region
  - GGGenome API parameter: `/1/` in URL

**Strand**:
- Both **plus (+) and minus (-) strands** are searched
- Plus strand: direct sequence match
- Minus strand: reverse complement match
- PAM recognition differs by strand (NGG on plus, CCN on minus)

### Annotation Parameters

**Annotation Database**:
- **GENCODE V39** (default)
- Exon annotations: `UCSC_exons_modif_canonical.bed`
- Intron annotations: `UCSC_introns_modif_canonical.bed`
- Files are UCSC-derived, modified BED format
- Can be customized via `-e` and `-i` options

**BED Format**:
- Standard 6-column BED format: `chrom`, `start`, `end`, `name`, `score`, `strand`
- Coordinates are 0-based, half-open (start inclusive, end exclusive)
- Used for genomic interval operations with bedtools

### GGGenome API Parameters

**API Base URL**: `https://gggenome.dbcls.jp/ja/{genome}/{mismatch}/{strand}/nogap/{sequence}.csv`

- `{genome}`: Genome assembly (e.g., hg38)
- `{mismatch}`: Maximum mismatches allowed (0-3)
- `{strand}`: Strand direction (+ or -)
- `nogap`: No gap alignment mode
- `{sequence}`: Query sequence (spacer+PAM or seed+PAM)

**CSV Format**:
- First 5 lines: Metadata/comments (skipped during parsing)
- Data columns: chrom, strand, start, end, snippet, snippet_pos, snippet_end, query, sbjct, align, edit, match, mis, del, ins
- Key columns:
  - Column 1: Chromosome
  - Column 2: Strand
  - Column 3: Start position
  - Column 4: End position
  - Column 9: Subject sequence (sbjct) - the matched genomic sequence

### Primer-BLAST Parameters

**Primer Design Window**:
- Upstream: 500bp before start, 150bp before end (5' primer region)
- Downstream: 150bp after start, 500bp after end (3' primer region)
- Product size: 400-800bp
- These parameters are optimized for PCR validation of off-target sites

**RefSeq ID Mapping**:
- Chromosome names (chr1, chr2, etc.) are converted to RefSeq IDs (NC_000001.11, etc.)
- Required for Primer-BLAST API compatibility
- Mapping table includes all autosomes (1-22) and sex chromosomes (X, Y)

### Default Values Summary

| Parameter | Default | Description |
|-----------|---------|-------------|
| Spacer length | 20 nt | Standard gRNA length |
| Seed length | 12 nt | Options: 8 to 12 (inclusive) |
| PAM | NGG | SpCas9 PAM sequence |
| Genome | hg38 | Human genome assembly |
| Full mismatch | 3 | Maximum mismatches for full sequence |
| Seed mismatch | 1 | Maximum mismatches for seed sequence |
| Exon DB | `data/UCSC_exons_modif_canonical.bed` | Exon annotations |
| Intron DB | `data/UCSC_introns_modif_canonical.bed` | Intron annotations |
| Output directory | `analysis/` | Results output location |

## Project Structure

```
OT_mapper/
├── bash/                    # Bash implementation
│   ├── scripts/
│   │   ├── utils.sh        # Shared utilities (logging, validation, helpers)
│   │   ├── OT_detector.sh  # Complete pipeline (detection + annotation + primer)
│   │   ├── OT_mapper.sh    # Gene annotation mapping (can be run standalone)
│   │   └── primer_generate.R  # Primer-BLAST URL generation
│   └── data/               # Annotation BED files
│
├── R/                      # R implementation
│   ├── scripts/
│   │   ├── OT_mapper.R     # Core OT detection functions
│   │   ├── OT_annotator.R  # Gene annotation functions
│   │   ├── analysis.R      # Complete pipeline (command-line script)
│   │   ├── primer_generate.R  # Primer-BLAST URL generation
│   │   ├── gggenome_functions.R  # GGGenome API helpers
│   │   ├── server.R        # Shiny app server
│   │   └── ui.R            # Shiny app UI
│   └── data/               # Annotation BED files
│
└── README.md               # This file
```

### Code Organization

**Bash Scripts:**
- **`utils.sh`**: Shared utility functions for logging, validation, file operations, and common tasks
- **`OT_detector.sh`**: Complete pipeline script that runs:
  - Off-target candidate detection
  - Gene annotation mapping (automatic, can be skipped with `--skip-annotation`)
  - Primer-BLAST link generation (automatic, can be skipped with `--skip-primer`)
  - Modular design with clear function separation
- **`OT_mapper.sh`**: Standalone script for annotation mapping (can be run independently)

**R Scripts:**
- **`analysis.R`**: Main command-line script with comprehensive error handling
- **`OT_mapper.R`**: Core functions for GGGenome queries and overlap detection
- **`OT_annotator.R`**: Gene annotation functions using bedtools
- **`gggenome_functions.R`**: Helper functions for GGGenome API interactions

**Design Principles:**
- **Modularity**: Functions are small, focused, and reusable
- **Error Handling**: Comprehensive validation and clear error messages
- **Logging**: Consistent, user-friendly progress and status messages
- **Documentation**: Clear comments and function descriptions
- **Maintainability**: Easy to read, modify, and extend

## Workflow

Both implementations follow a similar workflow:

1. **OT Candidate Detection** (`OT_detector.sh` / `OT_mapper.R`)
   - Accepts 20nt gRNA spacer sequence (without PAM)
   - Automatically appends `NGG` PAM
   - Queries GGGenome for full sequence (3 mismatches) and seed sequence (1 mismatch)
   - Filters candidates based on seed overlap and exact PAM matching
   - Outputs BED file with candidate coordinates

2. **Gene Annotation** (`OT_mapper.sh` / `OT_annotator.R`)
   - Maps OT candidates to GENCODE V39 annotations
   - Identifies exon/intron overlaps using `bedtools intersect`
   - Outputs annotated TSV file

3. **Primer Link Generation** (`primer_generate.R`)
   - Converts chromosome coordinates to RefSeq IDs
   - Generates Primer-BLAST URLs for each candidate
   - Outputs TSV with Primer-BLAST links

## Installation

### Dependencies

**Common dependencies** (both versions):
- `bedtools` (v2.30+)
- `wget` (for Bash version)
- Internet connection (for GGGenome API queries)

**Bash version additional**:
- `bash` (v4.0+)
- `R` with `optparse` package

**R version additional**:
- `R` (v4.0+)
- R packages: `optparse`, `httr`, `readr`, `dplyr`, `stringr`, `DT` (for Shiny app), `GenomicRanges` (optional, for advanced features)

### Setup

#### Using Conda (Recommended)

```bash
# Create conda environment
conda create -n OT -c bioconda r-optparse bedtools wget r-httr r-readr r-dplyr r-stringr
conda activate OT

# Install additional R packages if needed
Rscript -e "install.packages(c('DT', 'shiny', 'shinydashboard'), repos='https://cran.rstudio.com/')"
```

#### Using Homebrew (macOS)

```bash
brew install bedtools
# R packages can be installed via R console or Rscript
```

## Usage

### Bash Version

#### Quick Start (Interactive Mode)

```bash
cd bash
bash scripts/OT_detector.sh
```

The script will prompt you interactively:
- **Spacer sequence**: 20nt gRNA sequence (without PAM), e.g., `GCTGAAGCACTGCACGCCGT`
- **Seed length**: Choose `8` to `12` nucleotides (inclusive)

**Note**: By default, the script runs the complete pipeline:
1. Off-target candidate detection
2. Gene annotation mapping (automatic)
3. Primer-BLAST link generation (automatic)

Use `--skip-annotation` or `--skip-primer` to skip specific steps if needed.

#### Command-Line Mode

```bash
# Basic usage
bash scripts/OT_detector.sh -s GCTGAAGCACTGCACGCCGT -l 12

# With custom options
bash scripts/OT_detector.sh \
  -s GCTGAAGCACTGCACGCCGT \
  -l 12 \
  -p NGG \
  -g hg38 \
  -f 3 \
  -m 1 \
  -o results

# Show help
bash scripts/OT_detector.sh -h
```

**OT_detector.sh Options**:
- `-s, --spacer`: 20nt gRNA spacer sequence (without PAM) [required]
- `-l, --seed-length`: Seed sequence length: 8 to 12 (inclusive) [default: 12]
- `-p, --pam`: PAM sequence [default: NGG]
- `-g, --genome`: Genome assembly [default: hg38]
- `-f, --full-mismatch`: Mismatch tolerance for full sequence [default: 3]
- `-m, --seed-mismatch`: Mismatch tolerance for seed sequence [default: 1]
- `-o, --output-dir`: Output directory [default: analysis]
- `--skip-annotation`: Skip gene annotation step [default: false]
- `--skip-primer`: Skip Primer-BLAST link generation [default: false]
- `-h, --help`: Show help message

#### Output Files

- `{output_dir}/OT/OT_candidate.bed`: OT candidate BED file
- `{output_dir}/OT/OT_list_final.csv`: Full CSV data for candidates
- `{output_dir}/OT/UCSC_list_final.csv`: Sequence and coordinate list
- `{output_dir}/OT_mapper_results/OT_mapped.tsv`: Annotated TSV with exon/intron information
- `{output_dir}/OT_mapper_results/OT_with_primer.tsv`: Final TSV with Primer-BLAST URLs

#### Complete Pipeline (Recommended)

```bash
# Run complete pipeline (detection + annotation + primer generation)
bash scripts/OT_detector.sh -s GCTGAAGCACTGCACGCCGT -l 12

# Skip primer generation if not needed
bash scripts/OT_detector.sh -s GCTGAAGCACTGCACGCCGT -l 12 --skip-primer

# Skip both annotation and primer generation (detection only)
bash scripts/OT_detector.sh -s GCTGAAGCACTGCACGCCGT -l 12 --skip-annotation
```

#### Manual Step-by-Step Execution (Advanced)

If you need to run steps separately or with custom options:

```bash
# Step 1: Detect OT candidates only
bash scripts/OT_detector.sh -s GCTGAAGCACTGCACGCCGT -l 12 --skip-annotation

# Step 2: Map to gene annotations (standalone)
bash scripts/OT_mapper.sh analysis/OT/OT_candidate.bed

# Step 2: Map with custom annotation files
bash scripts/OT_mapper.sh analysis/OT/OT_candidate.bed \
  -e custom_exons.bed \
  -i custom_introns.bed \
  -o custom_output

# Step 3: Generate Primer-BLAST links (standalone)
Rscript scripts/primer_generate.R \
  -i analysis/OT_mapper_results/OT_mapped.tsv \
  -o analysis/OT_mapper_results/OT_with_primer.tsv
```

**OT_mapper.sh Options**:
- `-e, --exon-db`: Exon annotation BED file [default: data/UCSC_exons_modif_canonical.bed]
- `-i, --intron-db`: Intron annotation BED file [default: data/UCSC_introns_modif_canonical.bed]
- `-o, --output-dir`: Output directory [default: analysis/OT_mapper_results]
- `-h, --help`: Show help message

### R Version

#### Option 1: Interactive R/RStudio

```r
# Source the required scripts
source("scripts/OT_mapper.R")
source("scripts/OT_annotator.R")

# 1. GGGenome search (PAM=NGG)
list <- gggenome_to_dataframe("TCGCCCAGCGACCCTGCTCC", 8)

# 2. PAM exact match filtering
filtered_list <- filter_exact_pam(list)

# 3. Overlap detection
overlaps <- find_overlaps(filtered_list)

# 4. Combine results
combined_df <- combine_results(overlaps)

# 5. Gene annotation (exon/intron mapping)
annotated_df <- annotate_with_bedtools(
  combined_df,
  exon_db = "data/UCSC_exons_modif_canonical.bed",
  intron_db = "data/UCSC_introns_modif_canonical.bed",
  output_file = "annotated_offtargets.tsv"
)
```

#### Option 2: Command Line (Complete Pipeline)

```bash
cd R

# Complete pipeline (detection + annotation + primer generation)
Rscript scripts/analysis.R -s TCGCCCAGCGACCCTGCTCC -l 12

# Skip primer generation if not needed
Rscript scripts/analysis.R -s TCGCCCAGCGACCCTGCTCC -l 12 --skip-primer

# With all options
Rscript scripts/analysis.R \
  -s TCGCCCAGCGACCCTGCTCC \
  -l 12 \
  -p NGG \
  -g hg38 \
  -f 3 \
  -m 1 \
  -e scripts/data/UCSC_exons_modif_canonical.bed \
  -i scripts/data/UCSC_introns_modif_canonical.bed \
  -o results.tsv

# Legacy positional arguments (still supported)
Rscript scripts/analysis.R TCGCCCAGCGACCCTGCTCC 12 NGG results.tsv

# Show help
Rscript scripts/analysis.R -h
```

**Note**: By default, `analysis.R` runs the complete pipeline:
1. Off-target candidate detection
2. Gene annotation mapping (automatic)
3. Primer-BLAST link generation (automatic)

**Command-Line Options**:
- `-s, --spacer`: 20nt gRNA spacer sequence (without PAM) [required]
- `-l, --seed-length`: Seed sequence length: 8 to 12 (inclusive) [default: 12]
- `-p, --pam`: PAM sequence [default: NGG]
- `-g, --genome`: Genome assembly [default: hg38]
- `-f, --full-mismatch`: Mismatch tolerance for full sequence [default: 3]
- `-m, --seed-mismatch`: Mismatch tolerance for seed sequence [default: 1]
- `-e, --exon-db`: Exon annotation BED file [default: scripts/data/UCSC_exons_modif_canonical.bed]
- `-i, --intron-db`: Intron annotation BED file [default: scripts/data/UCSC_introns_modif_canonical.bed]
- `-o, --output`: Output file path [default: annotated_offtargets.tsv]
- `--skip-primer`: Skip Primer-BLAST link generation [default: false]
- `-h, --help`: Show help message

**Output Files**:
- `{output_file}`: Annotated TSV with exon/intron information
- `{output_file}_with_primer.tsv`: Final TSV with Primer-BLAST URLs (if not skipped)

#### Option 3: Shiny Web Application

**Starting the Shiny App**:

**Method 1: Using app.R (Recommended)**:
```bash
cd R/scripts
Rscript -e "shiny::runApp(port=3838)"
```

**Method 2: Direct runApp**:
```bash
cd R/scripts
Rscript -e "shiny::runApp('.', port=3838)"
```

**Method 3: In R/RStudio**:
```r
# Set working directory to R/scripts
setwd("R/scripts")

# Load Shiny and run the app
library(shiny)
runApp(port = 3838)
# Or simply: runApp() if app.R exists
```

Then open your browser to `http://localhost:3838` (or the URL displayed in the console).

**Shiny App Features**:
- Interactive input form for spacer sequence, seed length, and PAM
- Real-time progress updates during analysis
- Results table with sortable columns
- Download options:
  - Summary table (CSV)
  - Full results (CSV)
  - Annotated results (TSV)
- Summary statistics display

**Requirements for Shiny App**:
- R packages: `shiny`, `shinydashboard`, `DT`, `httr`, `readr`, `dplyr`, `stringr`, `GenomicRanges`
- Internet connection (for GGGenome API)
- bedtools (for overlap detection and annotation)
- Annotation database files in `scripts/data/` directory

**Install Shiny packages** (if not already installed):
```r
install.packages(c("shiny", "shinydashboard", "DT"))
```

## Input/Output Details

### Input

**Interactive Input**:
- **gRNA spacer**: 20 nucleotides, without PAM (PAM is automatically appended as `NGG`)
- **Seed length**: 8 to 12 nucleotides (inclusive)

**Example**:
- Spacer: `GCTGAAGCACTGCACGCCGT`
- Seed: `12`

### Annotation Reference Data

The tool uses GENCODE V39 annotation files (UCSC-derived, modified BED format):
- `UCSC_exons_modif_canonical.bed`: Exon annotations
- `UCSC_introns_modif_canonical.bed`: Intron annotations

These files are stored in `bash/data/` and `R/scripts/data/` directories.

### Output Files

1. **OT_candidate.bed**: BED format with chromosome, start, end coordinates
2. **OT_mapped.tsv**: Tab-separated file with:
   - Chromosome, start, end coordinates
   - Gene annotation (exon/intron)
   - Additional metadata from bedtools intersect
3. **OT_with_primer.tsv**: Same as OT_mapped.tsv with Primer-BLAST URLs added

## Important Notes

### Network Requirement
- **Internet connection required**: OT detection queries GGGenome API via `wget` (Bash) or `httr` (R)
- GGGenome API: `https://gggenome.dbcls.jp/hg38/`

### PAM Limitation
- Currently supports **NGG PAM only** (SpCas9) in the Bash version
- R version has more flexible PAM handling (configurable PAM parameter)
- Both versions filter for exact PAM matches to ensure biological relevance

### Annotation Filtering
- **bedtools removes candidates without annotations**: Candidates that don't overlap with exons or introns will not appear in `OT_mapped.tsv`
- Always check `OT_candidate.bed` to see all detected candidates, including those without gene annotations

### Column Consistency
- `primer_generate.R` expects chromosome in the first column (e.g., `chr1`)
- The Bash version's `OT_mapper.sh` uses `cut -f 5,1-3,8-9`, which may affect column order
- If primer generation fails, verify that columns 1-3 of `OT_mapped.tsv` are `chr/start/end`

### Output Directory
- Results are saved in `analysis/` subdirectory
- Intermediate files are stored in `analysis/intermediate/`
- No timestamps are added to output filenames (results overwrite previous runs)

## Example Usage

### Example 1: Bash Version (Interactive)

```bash
cd bash
bash scripts/OT_detector.sh
# Input when prompted:
# Spacer: TCGCCCAGCGACCCTGCTCC
# Seed: 12
# The script will automatically run detection, annotation, and primer generation
```

### Example 1b: Bash Version (Command-Line)

```bash
cd bash
# Complete pipeline (automatic)
bash scripts/OT_detector.sh -s TCGCCCAGCGACCCTGCTCC -l 12
```

### Example 2: R Command Line

```bash
cd R
# Complete pipeline (automatic)
Rscript scripts/analysis.R -s TCGCCCAGCGACCCTGCTCC -l 12

# Legacy positional arguments (still supported)
Rscript scripts/analysis.R TCGCCCAGCGACCCTGCTCC 12 NGG my_results.tsv
```

### Example 3: R Shiny App

```bash
cd R/scripts
Rscript -e "shiny::runApp(port=3838)"
# Then open http://localhost:3838 in your browser
```

### Example 3: R Interactive

```r
source("scripts/OT_mapper.R")
source("scripts/OT_annotator.R")

# Run analysis
list <- gggenome_to_dataframe("GCTGAAGCACTGCACGCCGT", 12, "NGG")
filtered <- filter_exact_pam(list, "NGG")
overlaps <- find_overlaps(filtered)
combined <- combine_results(overlaps)
annotated <- annotate_with_bedtools(
  combined,
  exon_db = "data/UCSC_exons_modif_canonical.bed",
  intron_db = "data/UCSC_introns_modif_canonical.bed",
  output_file = "results.tsv"
)
```

## Version Comparison

| Feature | Bash Version | R Version |
|---------|-------------|-----------|
| Entry Point | `scripts/OT_detector.sh` | `analysis.R`, Shiny app |
| Interactive Input | Yes | Yes (R console) + Shiny app |
| Command Line | Limited | Full support |
| Overlap Detection | **bedtools intersect** (coordinate-based) | **bedtools intersect** (coordinate-based) |
| PAM Flexibility | NGG only | Configurable |
| Web Interface | No | Yes (Shiny) |
| Progress Tracking | Basic | Advanced (Shiny) |
| Error Handling | Basic | Advanced |

**Note**: Both versions now use the same coordinate-based overlap detection algorithm using `bedtools intersect`, ensuring consistent results between implementations.

## Troubleshooting

### Common Issues

1. **"bedtools not found"**
   - Install bedtools: `conda install bedtools` or `brew install bedtools`
   - Verify: `which bedtools`

2. **"GGGenome download failed"**
   - Check internet connection
   - Verify GGGenome API is accessible: `curl https://gggenome.dbcls.jp/hg38/3/+/nogap/TCGCCCAGCGACCCTGCTCCNGG.csv`

3. **"No annotations found"**
   - Verify annotation BED files exist in `data/` directory
   - Check file paths in scripts
   - Some candidates may not overlap with exons/introns (check `OT_candidate.bed`)

4. **"Primer generation error"**
   - Verify `OT_mapped.tsv` column order (chr, start, end in first 3 columns)
   - Check RefSeq ID conversion table in `primer_generate.R`

## Contributions and Citations

### Development

This tool was developed by **Ryo** and **Gabriel** for CRISPR off-target analysis.

### Off-Target Definition

The off-target detection methodology and definitions used in this software were established by:
- **Dr. Knut Woltjen**
- **Dr. Ryo Niwa**

### Acknowledgments

- **GGGenome** (https://gggenome.dbcls.jp/) for off-target search API
- **GENCODE** for gene annotations
- **UCSC** for genome browser and annotation resources

