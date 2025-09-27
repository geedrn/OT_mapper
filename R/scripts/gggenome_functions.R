# gggenome_functions.R

library(httr)
library(readr)
library(dplyr)
library(stringr)

# スペース除去とバリデーション用の関数を改善
clean_and_validate_input <- function(input_string) {
  # すべての空白文字を除去
  cleaned <- gsub("[[:space:]]", "", input_string)
  return(cleaned)
}

# Function to query GGGenome API
# GGGenomeクエリ関数の改善
query_gggenome <- function(sequence, strand, mismatch) {
  base_url <- "https://gggenome.dbcls.jp/hg38/"
  encoded_sequence <- URLencode(sequence, reserved = TRUE)
  url <- paste0(base_url, mismatch, "/", strand, "/", encoded_sequence, ".csv")
  
  print(paste("Querying URL:", url))
  
  response <- GET(url)
  stop_for_status(response)
  
  result <- content(response, "parsed")
  
  if (is.null(result) || length(result) == 0) {
    return(list())  # 空のリストを返す
  }
  
  return(result)
}

# Function to process GGGenome results
process_gggenome <- function(content, sequence) {
  if (is.null(content) || content == "") {
    return(data.frame())
  }
  
  # Split content into lines
  lines <- strsplit(content, "\n")[[1]]
  
  # Extract metadata
  metadata <- lines[1:5]
  
  # Extract column names
  col_names <- strsplit(str_remove(metadata[5], "# "), ",")[[1]]
  
  # Parse data rows
  data_rows <- lines[6:length(lines)]
  data <- read_csv(paste(data_rows, collapse = "\n"), 
                   col_names = col_names, 
                   show_col_types = FALSE)
  
  # Calculate mismatches if necessary columns exist
  if (all(c("mis", "del", "ins") %in% names(data))) {
    data <- data %>%
      mutate(mismatches = as.numeric(mis) + as.numeric(del) + as.numeric(ins))
  } else {
    data$mismatches <- 0  # Default value if columns are missing
  }
  
  return(data)
}

