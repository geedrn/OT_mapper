#' GGGenomeからデータを取得してデータフレームに変換（PAM対応版）
#'
#' @param spacer gRNAスペーサー配列（例："TCGCCCAGCGACCCTGCTCC"）
#' @param seed_length シード配列の長さ
#' @param pam PAM配列（例："NGG"）
#' @param full_mismatch 全長配列の最大ミスマッチ数（デフォルト: 3）
#' @param seed_mismatch シード配列の最大ミスマッチ数（デフォルト: 1）
#' @param verbose 詳細出力を表示するか（デフォルト: TRUE、ShinyアプリではFALSE推奨）
#' @return 全長データとシード配列データのリスト（+/-鎖）
gggenome_to_dataframe <- function(spacer, seed_length, pam = "NGG", full_mismatch = 3, seed_mismatch = 1, verbose = TRUE) {
  if (missing(spacer) || !is.character(spacer) || length(spacer) != 1) {
    stop("spacerは1つの文字列である必要があります")
  }
  
  if (missing(seed_length) || !is.numeric(seed_length) || seed_length <= 0) {
    stop("seed_lengthは正の数である必要があります")
  }
  
  if (!is.character(pam) || length(pam) != 1) {
    stop("pamは1つの文字列である必要があります")
  }
  
  if (!is.numeric(full_mismatch) || full_mismatch < 0 || full_mismatch > 3) {
    stop("full_mismatchは0から3の整数である必要があります")
  }
  
  if (!is.numeric(seed_mismatch) || seed_mismatch < 0 || seed_mismatch > 3) {
    stop("seed_mismatchは0から3の整数である必要があります")
  }
  
  # ミスマッチ数を整数に変換
  full_mismatch <- as.integer(full_mismatch)
  seed_mismatch <- as.integer(seed_mismatch)
  
  # スペーサーとPAMを結合
  sequence_with_pam <- paste0(spacer, pam)
  
  # 詳細出力（Shinyアプリでは無効化）
  if (verbose) {
    cat("GGGenome検索を実行中:\n")
    cat("スペーサー: ", spacer, "\n")
    cat("PAM: ", pam, "\n")
    cat("検索配列: ", sequence_with_pam, "\n")
    cat("シード長: ", seed_length, "\n")
    cat("全長ミスマッチ許容数: ", full_mismatch, "\n")
    cat("シードミスマッチ許容数: ", seed_mismatch, "\n")
  }
  
  # シード配列を抽出（スペーサーの3'末端側）
  seed <- substr(spacer, nchar(spacer) - seed_length + 1, nchar(spacer))
  seed_with_pam <- paste0(seed, pam)
  
  if (verbose) {
    cat("抽出したシード配列: ", seed, "\n")
    cat("PAM付きシード配列: ", seed_with_pam, "\n")
  }
  
  # 一時ファイルのパスを設定
  temp_dir <- tempdir()
  plus_full_file <- file.path(temp_dir, "plus_full.csv")
  minus_full_file <- file.path(temp_dir, "minus_full.csv")
  plus_seed_file <- file.path(temp_dir, "plus_seed.csv")
  minus_seed_file <- file.path(temp_dir, "minus_seed.csv")
  
  if (verbose) {
    cat("ダウンロード先ディレクトリ: ", temp_dir, "\n")
    cat("ファイルパス:\n")
    cat("- 全長 + 鎖: ", plus_full_file, "\n")
    cat("- 全長 - 鎖: ", minus_full_file, "\n")
    cat("- シード + 鎖: ", plus_seed_file, "\n")
    cat("- シード - 鎖: ", minus_seed_file, "\n")
  }
  
  # GGGenome APIのURL（ヒトゲノムhg38、ミスマッチ設定）
  # PAMを含めた検索を行う
  plus_full_url <- paste0("https://gggenome.dbcls.jp/hg38/", full_mismatch, "/+/nogap/", sequence_with_pam, ".csv")
  minus_full_url <- paste0("https://gggenome.dbcls.jp/hg38/", full_mismatch, "/-/nogap/", sequence_with_pam, ".csv")
  
  # シード配列のAPI URL（ミスマッチ許容数を指定、PAMを含む）
  plus_seed_url <- paste0("https://gggenome.dbcls.jp/hg38/", seed_mismatch, "/+/nogap/", seed_with_pam, ".csv")
  minus_seed_url <- paste0("https://gggenome.dbcls.jp/hg38/", seed_mismatch, "/-/nogap/", seed_with_pam, ".csv")
  
  if (verbose) {
    cat("GGGenome APIを呼び出し中...\n")
    cat("【全長配列検索】\n")
    cat("+ 鎖 URL: ", plus_full_url, "\n")
    cat("- 鎖 URL: ", minus_full_url, "\n")
  }
  
  # 全長データのダウンロード
  download_success_plus_full <- tryCatch({
    utils::download.file(plus_full_url, plus_full_file, quiet = !verbose)
    TRUE
  }, error = function(e) {
    if (verbose) cat("+ 鎖全長データのダウンロードに失敗しました: ", e$message, "\n")
    FALSE
  })
  
  download_success_minus_full <- tryCatch({
    utils::download.file(minus_full_url, minus_full_file, quiet = !verbose)
    TRUE
  }, error = function(e) {
    if (verbose) cat("- 鎖全長データのダウンロードに失敗しました: ", e$message, "\n")
    FALSE
  })
  
  if (verbose) {
    cat("【シード配列検索】\n")
    cat("+ 鎖 URL: ", plus_seed_url, "\n")
    cat("- 鎖 URL: ", minus_seed_url, "\n")
  }
  
  # シード配列データのダウンロード
  download_success_plus_seed <- tryCatch({
    utils::download.file(plus_seed_url, plus_seed_file, quiet = !verbose)
    TRUE
  }, error = function(e) {
    if (verbose) cat("+ 鎖シードデータのダウンロードに失敗しました: ", e$message, "\n")
    FALSE
  })
  
  download_success_minus_seed <- tryCatch({
    utils::download.file(minus_seed_url, minus_seed_file, quiet = !verbose)
    TRUE
  }, error = function(e) {
    if (verbose) cat("- 鎖シードデータのダウンロードに失敗しました: ", e$message, "\n")
    FALSE
  })
  
  # ダウンロードしたCSVファイルを処理する関数
  process_gggenome_file <- function(file_path, file_type = "", verbose_param = verbose) {
    if (!file.exists(file_path) || file.size(file_path) == 0) {
      if (verbose_param) {
        if (file_type != "") {
          cat("ファイルが存在しないか空です (", file_type, "): ", basename(file_path), "\n")
        } else {
          cat("ファイルが存在しないか空です: ", file_path, "\n")
        }
      }
      return(data.frame())
    }
    
    # ファイル全体を読み込む
    lines <- readLines(file_path)
    
    # コメント行を除外
    data_lines <- lines[!grepl("^#", lines)]
    
    if (length(data_lines) == 0) {
      if (verbose_param) {
        if (file_type != "") {
          cat("データ行が見つかりません (", file_type, "): ", basename(file_path), "\n")
        } else {
          cat("データ行が見つかりません: ", basename(file_path), "\n")
        }
      }
      return(data.frame())
    }
    
    # データフレームに変換
    result <- tryCatch({
      # タブをカンマに置き換えてcsvとして読み込む
      temp_file <- tempfile(fileext = ".csv")
      writeLines(gsub("\t", ",", data_lines), temp_file)
      df <- read.csv(temp_file, header = FALSE, stringsAsFactors = FALSE)
      unlink(temp_file)
      
      # 列名を設定
      col_names <- c("chrom", "strand", "start", "end", "snippet", "snippet_pos", 
                     "snippet_end", "query", "sbjct", "align", "edit", "match", 
                     "mis", "del", "ins")
      
      colnames(df) <- col_names[seq_len(min(ncol(df), length(col_names)))]
      
      # 数値型に変換
      for (col in c("start", "end", "match", "mis", "del", "ins")) {
        if (col %in% colnames(df)) {
          df[[col]] <- suppressWarnings(as.integer(df[[col]]))
        }
      }
      
      # 位置情報を追加
      if (all(c("chrom", "start", "end") %in% colnames(df))) {
        df$location <- paste0(df$chrom, ":", df$start, "-", df$end)
      }
      
      df
    }, error = function(e) {
      if (verbose_param) {
        if (file_type != "") {
          cat("ファイルの解析に失敗しました (", file_type, "): ", e$message, "\n")
        } else {
          cat("ファイルの解析に失敗しました: ", e$message, "\n")
        }
      }
      return(data.frame())
    })
    
    return(result)
  }
  
  # ダウンロードした結果を処理
  plus_full_df <- if (download_success_plus_full) process_gggenome_file(plus_full_file, "全長 + 鎖", verbose) else data.frame()
  minus_full_df <- if (download_success_minus_full) process_gggenome_file(minus_full_file, "全長 - 鎖", verbose) else data.frame()
  plus_seed_df <- if (download_success_plus_seed) process_gggenome_file(plus_seed_file, "シード + 鎖", verbose) else data.frame()
  minus_seed_df <- if (download_success_minus_seed) process_gggenome_file(minus_seed_file, "シード - 鎖", verbose) else data.frame()
  
  # 結果件数の表示
  if (verbose) {
    cat("【検索結果】\n")
    cat("全長配列 + 鎖の結果: ", nrow(plus_full_df), " 件\n")
    cat("全長配列 - 鎖の結果: ", nrow(minus_full_df), " 件\n")
    cat("シード配列 + 鎖の結果: ", nrow(plus_seed_df), " 件\n")
    cat("シード配列 - 鎖の結果: ", nrow(minus_seed_df), " 件\n")
  }
  
  # 結果をリストとして返す
  result <- list(
    plus_full = plus_full_df,
    minus_full = minus_full_df,
    plus_seed = plus_seed_df,
    minus_seed = minus_seed_df,
    spacer = spacer,
    seed = seed,
    pam = pam,
    full_sequence = sequence_with_pam,
    seed_sequence = seed_with_pam
  )
  
  return(result)
}

#' シード配列と全長配列のオーバーラップを抽出する関数（PAM対応版）
#'
#' @param list PAM対応版gggenome_to_dataframe関数の結果リスト
#' @param use_bedtools bedtoolsを使用するかどうか
#' @return オーバーラップ結果を含むリスト
#' シード配列と全長配列のオーバーラップを抽出する関数（PAM対応版、全長配列出力）
#'
#' @param list PAM対応版gggenome_to_dataframe関数の結果リスト
#' @param use_bedtools bedtoolsを使用するかどうか
#' @return オーバーラップ結果を含むリスト（全長配列データを出力）
find_overlaps <- function(list, use_bedtools = TRUE, verbose = TRUE) {
  # 入力チェック
  if (!is.logical(verbose) || length(verbose) != 1) {
    verbose <- TRUE
    warning("verboseは論理値である必要があります。デフォルト値TRUEを使用します。")
  }
  
  if (!is.logical(use_bedtools) || length(use_bedtools) != 1) {
    use_bedtools <- TRUE
    warning("use_bedtoolsは論理値である必要があります。デフォルト値TRUEを使用します。")
  }
  
  if (verbose) cat("シード配列と全長配列のオーバーラップを検出中...\n")
  
  # 入力チェック
  if (!is.list(list) || !all(c("plus_full", "minus_full", "plus_seed", "minus_seed") %in% names(list))) {
    stop("listはgggenome_to_dataframe関数の結果である必要があります")
  }
  
  # データフレームが存在することを確認
  if (is.null(list$plus_full) || !is.data.frame(list$plus_full)) {
    list$plus_full <- data.frame()
  }
  if (is.null(list$minus_full) || !is.data.frame(list$minus_full)) {
    list$minus_full <- data.frame()
  }
  if (is.null(list$plus_seed) || !is.data.frame(list$plus_seed)) {
    list$plus_seed <- data.frame()
  }
  if (is.null(list$minus_seed) || !is.data.frame(list$minus_seed)) {
    list$minus_seed <- data.frame()
  }
  
  # bedtoolsのインストール確認
  bedtools_installed <- use_bedtools && (system("which bedtools", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0)
  if (use_bedtools && !bedtools_installed) {
    if (verbose) cat("bedtoolsが見つかりません。Rの内部処理でオーバーラップを検出します。\n")
  }
  
  # 処理用の一時ディレクトリ
  temp_dir <- tempdir()
  
  # 1. プラス鎖のオーバーラップ検出
  if (verbose) cat("プラス鎖のオーバーラップを検出中...\n")
  plus_overlaps <- data.frame()
  
  # 安全にnrow()をチェック
  plus_full_nrow <- if (is.data.frame(list$plus_full)) nrow(list$plus_full) else 0
  plus_seed_nrow <- if (is.data.frame(list$plus_seed)) nrow(list$plus_seed) else 0
  
  if (plus_full_nrow > 0 && plus_seed_nrow > 0) {
    # 必要な列が存在することを確認
    required_cols_full <- c("chrom", "start", "end", "strand")
    required_cols_seed <- c("chrom", "start", "end", "strand")
    
    has_required_cols_full <- all(required_cols_full %in% colnames(list$plus_full))
    has_required_cols_seed <- all(required_cols_seed %in% colnames(list$plus_seed))
    
    if (!has_required_cols_full) {
      if (verbose) cat("警告: プラス鎖全長データに必要な列がありません\n")
    } else if (!has_required_cols_seed) {
      if (verbose) cat("警告: プラス鎖シードデータに必要な列がありません\n")
    } else if (bedtools_installed) {
      # BEDファイル作成
      plus_full_bed <- file.path(temp_dir, "plus_full.bed")
      plus_seed_bed <- file.path(temp_dir, "plus_seed.bed")
      plus_overlap_bed <- file.path(temp_dir, "plus_overlap.bed")
      
      # 全長配列をBED形式で保存
      write.table(
        data.frame(
          list$plus_full$chrom,
          list$plus_full$start,
          list$plus_full$end,
          paste0("full_", 1:plus_full_nrow),
          rep(0, plus_full_nrow),
          list$plus_full$strand
        ),
        plus_full_bed,
        quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t"
      )
      
      # シード配列をBED形式で保存
      write.table(
        data.frame(
          list$plus_seed$chrom,
          list$plus_seed$start,
          list$plus_seed$end,
          paste0("seed_", 1:plus_seed_nrow),
          rep(0, plus_seed_nrow),
          list$plus_seed$strand
        ),
        plus_seed_bed,
        quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t"
      )
      
      # bedtools intersectで重なりを検出
      # -bを先に指定して-waを使用することで、全長配列(plus_full_bed)を出力
      cmd <- paste("bedtools intersect -b", plus_seed_bed, "-a", plus_full_bed, "-wa >", plus_overlap_bed)
      system(cmd)
      
      # 結果読み込み
      if (file.exists(plus_overlap_bed) && file.size(plus_overlap_bed) > 0) {
        overlap_bed <- read.table(plus_overlap_bed, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
        colnames(overlap_bed) <- c("chrom", "start", "end", "name", "score", "strand")
        
        # 元のデータフレームから一致する行を抽出
        idx <- numeric(0)
        for (i in seq_len(nrow(overlap_bed))) {
          matches <- which(
            list$plus_full$chrom == overlap_bed$chrom[i] &
              list$plus_full$start == overlap_bed$start[i] &
              list$plus_full$end == overlap_bed$end[i] &
              list$plus_full$strand == overlap_bed$strand[i]
          )
          idx <- c(idx, matches)
        }
        
        if (length(idx) > 0) {
          plus_overlaps <- list$plus_full[idx, ]
        }
      }
    } else {
      # Rでの処理: メモリを節約するために分割処理
      batch_size <- 1000
      batches <- ceiling(nrow(list$plus_seed) / batch_size)
      
      for (b in 1:batches) {
        start_idx <- (b - 1) * batch_size + 1
        end_idx <- min(b * batch_size, nrow(list$plus_seed))
        
        # 現在のバッチ
        current_batch <- list$plus_seed[start_idx:end_idx, ]
        
        # オーバーラップ検出 - 全長配列を収集
        for (i in seq_len(nrow(current_batch))) {
          # 重なりのある全長配列のインデックスを見つける
          overlap_idx <- which(
            list$plus_full$chrom == current_batch$chrom[i] &
              (
                (list$plus_full$start <= current_batch$start[i] & list$plus_full$end >= current_batch$start[i]) |
                  (list$plus_full$start <= current_batch$end[i] & list$plus_full$end >= current_batch$end[i]) |
                  (list$plus_full$start >= current_batch$start[i] & list$plus_full$end <= current_batch$end[i])
              ) &
              list$plus_full$strand == current_batch$strand[i]
          )
          
          if (length(overlap_idx) > 0) {
            # 重なりのある全長配列を追加
            plus_overlaps <- rbind(plus_overlaps, list$plus_full[overlap_idx, ])
          }
        }
        
        # 進捗表示
        if (verbose && (b %% 10 == 0 || b == batches)) {
          cat(sprintf("  プラス鎖 - バッチ処理: %d/%d 完了\n", b, batches))
        }
      }
      
      # 重複を削除（同じ全長配列が複数のシード配列と重なる場合）
      if (nrow(plus_overlaps) > 0) {
        plus_overlaps <- unique(plus_overlaps)
      }
    }
  }
  
  # 2. マイナス鎖のオーバーラップ検出
  if (verbose) cat("マイナス鎖のオーバーラップを検出中...\n")
  minus_overlaps <- data.frame()
  
  # 安全にnrow()をチェック
  minus_full_nrow <- if (is.data.frame(list$minus_full)) nrow(list$minus_full) else 0
  minus_seed_nrow <- if (is.data.frame(list$minus_seed)) nrow(list$minus_seed) else 0
  
  if (minus_full_nrow > 0 && minus_seed_nrow > 0) {
    # 必要な列が存在することを確認
    required_cols_full <- c("chrom", "start", "end", "strand")
    required_cols_seed <- c("chrom", "start", "end", "strand")
    
    has_required_cols_full <- all(required_cols_full %in% colnames(list$minus_full))
    has_required_cols_seed <- all(required_cols_seed %in% colnames(list$minus_seed))
    
    if (!has_required_cols_full) {
      if (verbose) cat("警告: マイナス鎖全長データに必要な列がありません\n")
    } else if (!has_required_cols_seed) {
      if (verbose) cat("警告: マイナス鎖シードデータに必要な列がありません\n")
    } else if (bedtools_installed) {
      # BEDファイル作成
      minus_full_bed <- file.path(temp_dir, "minus_full.bed")
      minus_seed_bed <- file.path(temp_dir, "minus_seed.bed")
      minus_overlap_bed <- file.path(temp_dir, "minus_overlap.bed")
      
      # 全長配列をBED形式で保存
      write.table(
        data.frame(
          list$minus_full$chrom,
          list$minus_full$start,
          list$minus_full$end,
          paste0("full_", 1:minus_full_nrow),
          rep(0, minus_full_nrow),
          list$minus_full$strand
        ),
        minus_full_bed,
        quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t"
      )
      
      # シード配列をBED形式で保存
      write.table(
        data.frame(
          list$minus_seed$chrom,
          list$minus_seed$start,
          list$minus_seed$end,
          paste0("seed_", 1:minus_seed_nrow),
          rep(0, minus_seed_nrow),
          list$minus_seed$strand
        ),
        minus_seed_bed,
        quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t"
      )
      
      # bedtools intersectで重なりを検出
      # -bを先に指定して-waを使用することで、全長配列(minus_full_bed)を出力
      cmd <- paste("bedtools intersect -b", minus_seed_bed, "-a", minus_full_bed, "-wa >", minus_overlap_bed)
      system(cmd)
      
      # 結果読み込み
      if (file.exists(minus_overlap_bed) && file.size(minus_overlap_bed) > 0) {
        overlap_bed <- read.table(minus_overlap_bed, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
        colnames(overlap_bed) <- c("chrom", "start", "end", "name", "score", "strand")
        
        # 元のデータフレームから一致する行を抽出
        idx <- numeric(0)
        for (i in seq_len(nrow(overlap_bed))) {
          matches <- which(
            list$minus_full$chrom == overlap_bed$chrom[i] &
              list$minus_full$start == overlap_bed$start[i] &
              list$minus_full$end == overlap_bed$end[i] &
              list$minus_full$strand == overlap_bed$strand[i]
          )
          idx <- c(idx, matches)
        }
        
        if (length(idx) > 0) {
          minus_overlaps <- list$minus_full[idx, ]
        }
      }
    } else {
      # Rでの処理: メモリを節約するために分割処理
      batch_size <- 1000
      batches <- ceiling(nrow(list$minus_seed) / batch_size)
      
      for (b in 1:batches) {
        start_idx <- (b - 1) * batch_size + 1
        end_idx <- min(b * batch_size, nrow(list$minus_seed))
        
        # 現在のバッチ
        current_batch <- list$minus_seed[start_idx:end_idx, ]
        
        # オーバーラップ検出 - 全長配列を収集
        for (i in seq_len(nrow(current_batch))) {
          # 重なりのある全長配列のインデックスを見つける
          overlap_idx <- which(
            list$minus_full$chrom == current_batch$chrom[i] &
              (
                (list$minus_full$start <= current_batch$start[i] & list$minus_full$end >= current_batch$start[i]) |
                  (list$minus_full$start <= current_batch$end[i] & list$minus_full$end >= current_batch$end[i]) |
                  (list$minus_full$start >= current_batch$start[i] & list$minus_full$end <= current_batch$end[i])
              ) &
              list$minus_full$strand == current_batch$strand[i]
          )
          
          if (length(overlap_idx) > 0) {
            # 重なりのある全長配列を追加
            minus_overlaps <- rbind(minus_overlaps, list$minus_full[overlap_idx, ])
          }
        }
        
        # 進捗表示
        if (verbose && (b %% 10 == 0 || b == batches)) {
          cat(sprintf("  マイナス鎖 - バッチ処理: %d/%d 完了\n", b, batches))
        }
      }
      
      # 重複を削除（同じ全長配列が複数のシード配列と重なる場合）
      if (nrow(minus_overlaps) > 0) {
        minus_overlaps <- unique(minus_overlaps)
      }
    }
  }
  
  # 結果表示
  if (verbose) {
    cat(sprintf("プラス鎖オーバーラップ: %d件（全長配列）\n", nrow(plus_overlaps)))
    cat(sprintf("マイナス鎖オーバーラップ: %d件（全長配列）\n", nrow(minus_overlaps)))
  }
  
  # 結果を返す
  result <- list(
    plus_overlaps = plus_overlaps,
    minus_overlaps = minus_overlaps,
    spacer = list$spacer,
    seed = list$seed,
    pam = list$pam
  )
  
  # 元のリストからメタデータがあれば追加
  metadata_fields <- c("full_sequence", "seed_sequence")
  for (field in metadata_fields) {
    if (!is.null(list[[field]])) {
      result[[field]] <- list[[field]]
    }
  }
  
  return(result)
}

#' PAMが完全一致する配列のみをフィルタリングする関数
#'
#' @param results gggenome_to_dataframe関数またはfind_overlaps関数の結果
#' @param pam 期待するPAM配列（デフォルト: "NGG"）
#' @param verbose 詳細出力を表示するか（デフォルト: TRUE）
#' @return PAMが完全一致する配列のみを含むリスト
filter_exact_pam <- function(results, pam = "NGG", verbose = TRUE) {
  if (verbose) {
    cat("PAMが完全一致する配列をフィルタリングしています...\n")
    cat("対象PAM: ", pam, "\n")
  }
  
  # 入力チェック
  if (!is.list(results)) {
    stop("resultsはリスト形式である必要があります")
  }
  
  # 入力タイプを判別
  is_overlaps_result <- all(c("plus_overlaps", "minus_overlaps") %in% names(results))
  is_gggenome_result <- all(c("plus_full", "minus_full", "plus_seed", "minus_seed") %in% names(results))
  
  if (!is_overlaps_result && !is_gggenome_result) {
    stop("サポートされていない結果形式です")
  }
  
  # PAMの長さを取得
  pam_length <- nchar(pam)
  
  # PAM位置のミスマッチをチェックする関数
  check_pam_match <- function(df, strand, pam, spacer_length = NULL, verbose_param = verbose) {
    if (nrow(df) == 0 || !("sbjct" %in% colnames(df))) {
      if (verbose_param) cat(strand, "鎖のデータが空または必要な列がありません\n")
      return(data.frame())
    }
    
    # sbjct列から直接PAM部分を抽出してチェック
    selected <- logical(nrow(df))
    
    # PAMパターンを正規表現に変換（Nを任意の塩基に）
    pam_pattern <- gsub("N", "[ATGC]", pam)
    
    # マイナス鎖の場合、PAMパターンも逆相補にする必要がある（一度だけ計算）
    pam_pattern_revcomp <- NULL
    if (strand == "-") {
      # NGG -> CCN (逆相補: 逆順にして相補塩基に変換)
      # 例: NGG -> GGN -> CCN
      pam_chars <- strsplit(pam, "")[[1]]
      pam_revcomp_chars <- character(length(pam_chars))
      for (k in seq_along(pam_chars)) {
        if (pam_chars[k] == "N") {
          pam_revcomp_chars[length(pam_chars) - k + 1] <- "N"
        } else {
          pam_revcomp_chars[length(pam_chars) - k + 1] <- switch(
            pam_chars[k],
            "A" = "T",
            "T" = "A",
            "G" = "C",
            "C" = "G",
            pam_chars[k]
          )
        }
      }
      pam_revcomp <- paste0(pam_revcomp_chars, collapse = "")
      pam_pattern_revcomp <- gsub("N", "[ATGC]", pam_revcomp)
    }
    
    for (i in seq_len(nrow(df))) {
      sbjct_seq <- toupper(df$sbjct[i])  # 大文字に変換
      sbjct_len <- nchar(sbjct_seq)
      
      if (sbjct_len < pam_length) {
        selected[i] <- FALSE
        next
      }
      
      # PAM部分を抽出
      if (strand == "+") {
        # プラス鎖: 配列末尾にPAM
        # spacer_lengthが指定されている場合はそれを使用、なければsbjctの末尾から
        if (!is.null(spacer_length) && spacer_length > 0) {
          if (sbjct_len >= spacer_length + pam_length) {
            pam_part <- substr(sbjct_seq, spacer_length + 1, spacer_length + pam_length)
          } else {
            selected[i] <- FALSE
            next
          }
        } else {
          # sbjctの末尾からPAM部分を抽出
          pam_part <- substr(sbjct_seq, sbjct_len - pam_length + 1, sbjct_len)
        }
      } else {
        # マイナス鎖: 配列先頭にPAM（逆相補配列なので、実際のPAMは逆相補になる）
        # マイナス鎖の場合、sbjctは逆相補配列なので、先頭がPAM部分
        pam_part <- substr(sbjct_seq, 1, pam_length)
      }
      
      # PAM部分がパターンに一致するかチェック
      if (strand == "+") {
        # プラス鎖: 正規表現でマッチング
        is_match <- grepl(pam_pattern, pam_part)
      } else {
        # マイナス鎖: 逆相補パターンでマッチング
        is_match <- grepl(pam_pattern_revcomp, pam_part)
      }
      
      selected[i] <- is_match
    }
    
    # フィルタリング結果の確認
    if (verbose_param) cat(strand, "鎖: ", sum(selected), "/", nrow(df), " 件が完全PAM一致\n")
    
    return(df[selected, ])
  }
  
  # スペーサー長を取得（全長配列のフィルタリング用）
  spacer_length <- NULL
  if (!is.null(results$spacer)) {
    spacer_length <- nchar(results$spacer)
  }
  
  # 結果構造に応じたフィルタリング
  if (is_overlaps_result) {
    # find_overlaps結果の処理
    plus_filtered <- check_pam_match(results$plus_overlaps, "+", pam, spacer_length = spacer_length, verbose_param = verbose)
    minus_filtered <- check_pam_match(results$minus_overlaps, "-", pam, spacer_length = spacer_length, verbose_param = verbose)
    
    # 元のメタデータを保持しつつ、フィルタリング結果を返す
    filtered_results <- results
    filtered_results$plus_overlaps <- plus_filtered
    filtered_results$minus_overlaps <- minus_filtered
    filtered_results$pam_filtered <- TRUE
  } else {
    # gggenome_to_dataframe結果の処理
    plus_full_filtered <- check_pam_match(results$plus_full, "+", pam, spacer_length = spacer_length, verbose_param = verbose)
    minus_full_filtered <- check_pam_match(results$minus_full, "-", pam, spacer_length = spacer_length, verbose_param = verbose)
    # シード配列の場合はspacer_lengthをNULLにして、末尾からPAMを抽出
    plus_seed_filtered <- check_pam_match(results$plus_seed, "+", pam, spacer_length = NULL, verbose_param = verbose)
    minus_seed_filtered <- check_pam_match(results$minus_seed, "-", pam, spacer_length = NULL, verbose_param = verbose)
    
    # 元のメタデータを保持しつつ、フィルタリング結果を返す
    filtered_results <- results
    filtered_results$plus_full <- plus_full_filtered
    filtered_results$minus_full <- minus_full_filtered
    filtered_results$plus_seed <- plus_seed_filtered
    filtered_results$minus_seed <- minus_seed_filtered
    filtered_results$pam_filtered <- TRUE
  }
  
  # 要約を表示
  if (verbose) cat("PAM完全一致フィルタリング完了\n")
  
  return(filtered_results)
}

#' CRISPR解析結果を1つのデータフレームに結合する
#'
#' @param results find_overlaps関数またはfilter_exact_pam関数の結果
#' @param include_metadata メタデータ列を含めるかどうか
#' @param verbose 詳細出力を表示するか（デフォルト: TRUE）
#' @return 結合されたデータフレーム
combine_results <- function(results, include_metadata = TRUE, verbose = TRUE) {
  if (verbose) cat("解析結果を1つのデータフレームに結合しています...\n")
  
  # 入力チェック
  if (!is.list(results)) {
    stop("resultsはリスト形式である必要があります")
  }
  
  # 結果タイプを判別
  is_overlaps_result <- all(c("plus_overlaps", "minus_overlaps") %in% names(results))
  
  if (!is_overlaps_result) {
    stop("サポートされていない結果形式です。find_overlapsまたはfilter_exact_pamの結果である必要があります")
  }
  
  # プラス鎖とマイナス鎖の結果を取得
  plus_data <- results$plus_overlaps
  minus_data <- results$minus_overlaps
  
  # 各データフレームにソース情報を追加
  if (nrow(plus_data) > 0) {
    plus_data$source <- "plus_strand"
  }
  
  if (nrow(minus_data) > 0) {
    minus_data$source <- "minus_strand"
  }
  
  # データフレームを結合
  combined_df <- rbind(plus_data, minus_data)
  
  # 結果が空の場合
  if (nrow(combined_df) == 0) {
    if (verbose) cat("結合されたデータフレームは空です\n")
    # 空のデータフレームを返す
    if (include_metadata) {
      return(data.frame(
        spacer = character(0),
        pam = character(0),
        result_count = integer(0),
        stringsAsFactors = FALSE
      ))
    } else {
      return(data.frame())
    }
  }
  
  # ミスマッチ数でソート（可能な場合）
  if ("mis" %in% colnames(combined_df)) {
    combined_df <- combined_df[order(combined_df$mis), ]
  }
  
  # メタデータを追加（必要な場合）
  if (include_metadata) {
    # メタデータ列を追加
    spacer <- if (!is.null(results$spacer)) results$spacer else "unknown"
    pam <- if (!is.null(results$pam)) results$pam else "unknown"
    
    # 行数を取得
    n <- nrow(combined_df)
    
    # メタデータ列を追加
    combined_df$spacer <- rep(spacer, n)
    combined_df$pam <- rep(pam, n)
    combined_df$result_id <- 1:n
  }
  
  # シーケンスの長さを追加（可能な場合）
  if ("sbjct" %in% colnames(combined_df)) {
    combined_df$sequence_length <- nchar(combined_df$sbjct)
  }
  
  if (verbose) cat(sprintf("結合完了: 合計 %d 件の結果\n", nrow(combined_df)))
  
  return(combined_df)
}