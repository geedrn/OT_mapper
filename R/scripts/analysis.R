#!/usr/bin/env Rscript

# CRISPRオフターゲット解析を実行するスクリプト
#
# 使用方法:
# Rscript run_crispr_analysis.R [配列] [シード長] [PAM] [出力ファイル]
#
# 例:
# Rscript run_crispr_analysis.R TCGCCCAGCGACCCTGCTCC 8 NGG results.tsv
#

# コマンドライン引数を解析
args <- commandArgs(trailingOnly = TRUE)

# デフォルト値の設定
spacer <- if (length(args) >= 1) args[1] else "TCGCCCAGCGACCCTGCTCC"
seed_length <- if (length(args) >= 2) as.numeric(args[2]) else 8
pam <- if (length(args) >= 3) args[3] else "NGG"
output_file <- if (length(args) >= 4) args[4] else "annotated_offtargets.tsv"

# 設定情報の表示
cat("=== CRISPR オフターゲット解析 ===\n")
cat("スペーサー配列: ", spacer, "\n")
cat("シード長: ", seed_length, "\n")
cat("PAM: ", pam, "\n")
cat("出力ファイル: ", output_file, "\n")
cat("\n")

# OT_mapper.RとOT_annotator.Rをソース
script_dir <- dirname(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))][1])
script_dir <- sub("--file=", "", script_dir)

if (script_dir == "") {
  script_dir <- getwd()
}

mapper_path <- file.path(script_dir, "OT_mapper.R")
annotator_path <- file.path(script_dir, "OT_annotator.R")

cat("スクリプトディレクトリ: ", script_dir, "\n")
cat("OT_mapper.R: ", mapper_path, "\n")
cat("OT_annotator.R: ", annotator_path, "\n")

source(mapper_path)
source(annotator_path)

# エクソン・イントロンデータベースのパス設定
# スクリプトディレクトリのサブディレクトリ "data" にあると仮定
exon_db <- file.path(script_dir, "data/UCSC_exons_modif_canonical.bed")
intron_db <- file.path(script_dir, "data/UCSC_introns_modif_canonical.bed")

cat("エクソンDB: ", exon_db, "\n")
cat("イントロンDB: ", intron_db, "\n")
cat("\n")

# 1. GGGenome検索（指定されたPAMを使用）
cat("ステップ1: GGGenomeでオフターゲット候補検索中...\n")
list <- gggenome_to_dataframe(spacer, seed_length, pam)

# 2. PAM完全一致フィルタリング
cat("\nステップ2: PAM完全一致でフィルタリング中...\n")
filtered_list <- filter_exact_pam(list, pam)

# 3. オーバーラップ検出
cat("\nステップ3: オーバーラップ検出中...\n")
overlaps <- find_overlaps(filtered_list)

# 4. 結果を1つのデータフレームに結合
cat("\nステップ4: 結果を結合中...\n")
combined_df <- combine_results(overlaps)

# 結果の確認
cat("\n検出されたオフターゲット候補: ", nrow(combined_df), " 件\n")

# 5. 遺伝子アノテーション（エクソン・イントロンマッピング）
cat("\nステップ5: 遺伝子アノテーション中...\n")
annotated_df <- annotate_with_bedtools(
  combined_df,
  exon_db = exon_db,
  intron_db = intron_db,
  output_file = output_file
)

# 最終結果の表示
cat("\n=== 解析完了 ===\n")
cat("オフターゲット解析結果: ", output_file, "\n")
cat("検出されたオフターゲット候補: ", nrow(combined_df), " 件\n")

if (is.null(annotated_df) || nrow(annotated_df) == 0) {
  cat("遺伝子アノテーションなし\n")
} else {
  cat("アノテーション行数: ", nrow(annotated_df), " 件\n")
}

cat("\n処理完了\n")