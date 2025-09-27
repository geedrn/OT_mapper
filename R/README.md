# OT mapper
## Introduction
This tool was developed by Ryo and Gabriel for identifying the OT candidates and annotating the gene information to the candidates. This tool can be controlled by the simple command line.

## Human annotation files
- GENCODE V39 annotation BED files were stored in `data/`.

## How to use
- On Rstudio

```R=
# OT_mapper.RとOT_annotator.Rをソース
source("OT_mapper.R")
source("OT_annotator.R")

# 1. GGGenome検索（PAM=NGGを使用）
list <- gggenome_to_dataframe("TCGCCCAGCGACCCTGCTCC", 8)

# 2. PAM完全一致フィルタリング
filtered_list <- filter_exact_pam(list)

# 3. オーバーラップ検出
overlaps <- find_overlaps(filtered_list)

# 4. 結果を1つのデータフレームに結合
combined_df <- combine_results(overlaps)

# 5. 結果の確認
print(paste("検出されたオフターゲット候補:", nrow(combined_df), "件"))
head(combined_df)

# 6. 遺伝子アノテーション（エクソン・イントロンマッピング）
annotated_df <- annotate_with_bedtools(
  combined_df,
  exon_db = "~/Documents/GitHub/OT_mapper/data/UCSC_exons_modif_canonical.bed",
  intron_db = "~/Documents/GitHub/OT_mapper/data/UCSC_introns_modif_canonical.bed",
  output_file = "annotated_offtargets.tsv"
)

# 7. アノテーション結果の確認
head(annotated_df)
```

- On terminal

```bash=
Rscript scripts/analysis.R TCGCCCAGCGACCCTGCTCC 8 NGG results.tsv
```

## Dependencies
```bash
# Look at homebrew website if you do not install them yet
# https://brew.sh/ja/
brew install bedtools
```

## Caution
Bedtools remove the candidates without annotatons. You need to make sure to check both combined_df and annotated_df to clarify all the candidates.
