#' オフターゲット候補を遺伝子アノテーションにマッピングする関数
#' 
#' オフターゲット候補を遺伝子アノテーションにマッピングする関数（修正版）
#' 
#' @param combined_df オフターゲット候補のデータフレーム
#' @param exon_db エクソンアノテーションBEDファイル
#' @param intron_db イントロンアノテーションBEDファイル
#' @param output_file 出力ファイルパス（デフォルトはNULL=出力なし）
#' @return アノテーション結果のデータフレーム
#' 
#' @details 出力データフレームの列構造:
#'   - V1-V6: 入力BEDファイルの列 (chrom, start, end, name, score, strand)
#'   - V7-V9: アノテーションBEDファイルの座標 (chrom, start, end)
#'   - V10: 遺伝子シンボル (geneSymbol)
#'   - V11: 特徴タイプ ("intron" または "exon")
#'   - V12: イントロン/エクソン番号 (0から始まる番号。例: 5 = 5番目のイントロン)
#'
annotate_with_bedtools <- function(combined_df, exon_db, intron_db, output_file = NULL) {
  # bedtoolsのインストール確認
  if (system("which bedtools", ignore.stdout = TRUE) != 0) {
    stop("bedtoolsがインストールされていないか、PATHに含まれていません")
  }
  
  # 必要な列の確認
  required_cols <- c("chrom", "start", "end")
  if (!all(required_cols %in% colnames(combined_df))) {
    stop("データフレームには chrom, start, end 列が必要です")
  }
  
  # ファイルの存在確認
  if (!file.exists(exon_db)) {
    stop(paste("エクソンDBファイルが見つかりません:", exon_db))
  }
  
  if (!file.exists(intron_db)) {
    stop(paste("イントロンDBファイルが見つかりません:", intron_db))
  }
  
  # 一時ファイルを作成
  temp_dir <- tempdir()
  temp_bed <- file.path(temp_dir, "offtargets.bed")
  temp_combined_db <- file.path(temp_dir, "combined_db.bed")
  temp_result <- file.path(temp_dir, "annotated_result.txt")
  
  # データフレームをBEDファイルに変換
  # 最低限の列を使用
  bed_data <- combined_df[, c("chrom", "start", "end")]
  
  # name, score, strand列の追加 (存在する場合は使用、なければ作成)
  if (!"name" %in% colnames(combined_df)) {
    bed_data$name <- paste0("OT_", 1:nrow(bed_data))
  } else {
    bed_data$name <- combined_df$name
  }
  
  if (!"score" %in% colnames(combined_df)) {
    bed_data$score <- 0
  } else {
    bed_data$score <- combined_df$score
  }
  
  if (!"strand" %in% colnames(combined_df)) {
    bed_data$strand <- "+"
  } else {
    bed_data$strand <- combined_df$strand
  }
  
  # BEDファイルに書き出し
  write.table(bed_data, temp_bed, 
              sep = "\t", 
              row.names = FALSE, 
              col.names = FALSE, 
              quote = FALSE)
  
  # エクソンとイントロンのDBを結合
  system(paste("cat", exon_db, intron_db, ">", temp_combined_db))
  
  # bedtools intersectを実行
  cmd <- paste("bedtools intersect -a", temp_bed, "-b", temp_combined_db, "-wa -wb >", temp_result)
  # Use print() instead of message() for Shiny compatibility
  if (!is.null(output_file)) {
    print(paste("実行コマンド: ", cmd))
  }
  
  status <- system(cmd)
  
  if (status != 0) {
    if (!is.null(output_file)) {
      print("Warning: bedtools実行中にエラーが発生しました")
    }
    return(NULL)
  }
  
  # 結果を読み込む
  if (file.exists(temp_result) && file.size(temp_result) > 0) {
    result_data <- read.table(temp_result, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    if (!is.null(output_file)) {
      print(paste("アノテーション完了:", nrow(result_data), "件"))
    }
    
    # 結果を出力（指定がある場合）
    if (!is.null(output_file)) {
      write.table(result_data, output_file, 
                  quote = FALSE, 
                  sep = "\t", 
                  row.names = FALSE, 
                  col.names = FALSE)
      print(paste("結果を保存しました:", output_file))
    }
    
    # 一時ファイルを削除
    unlink(c(temp_bed, temp_combined_db, temp_result))
    
    return(result_data)
  } else {
    if (!is.null(output_file)) {
      print("アノテーション結果は空でした")
    }
    unlink(c(temp_bed, temp_combined_db, temp_result))
    return(data.frame())
  }
}