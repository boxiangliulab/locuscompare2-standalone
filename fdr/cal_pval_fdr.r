# input file must contain "gene_id" column
# results are tab seperated table with 2 column (1 row): fdr, pvalue
library(qvalue)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("At least 2 arguments are rquired: input_file_path, output_file_path")
}

if (!require(qvalue)) {
  stop("qvalue not installed")
}

input_file_path = args[1]
output_file_path = args[2]
pval_col_name = args[3]
fdr_threshold = as.numeric(args[4])
if (is.na(input_file_path) || !file.exists(input_file_path) || !file.info(input_file_path)$size > 0) {
  print("input_file_path file does not exist or is empty!")
  q(save = "no")
}
if (is.na(pval_col_name) || tolower(pval_col_name) == 'na' || tolower(pval_col_name) == 'none') {
  pval_col_name = "pvalue"
}
if (is.na(fdr_threshold) || tolower(fdr_threshold) == 'na' || tolower(fdr_threshold) == 'none') {
  fdr_threshold = 0.05
}

data = read.table(file = input_file_path, header = TRUE, fill = TRUE)

# Filter out NA, -Inf, Inf, and out-of-range values from the p-value column
print("Filter1")
# 将 -Inf 替换为 0
data[[pval_col_name]][is.infinite(data[[pval_col_name]]) & data[[pval_col_name]] == -Inf] <- 0

print("Filter2")
# 过滤掉 NA 和非有限值，同时确保 p-value 在 0 到 1 之间
data <- data %>%
  filter(data[[pval_col_name]] >= 0 & data[[pval_col_name]] <= 1)

print(summary(data[[pval_col_name]]))

# 检查是否还有 NA 或 Inf
print(sum(is.na(data[[pval_col_name]])))  # 检查 NA 值的数量
print(sum(!is.finite(data[[pval_col_name]])))  # 检查 Inf 和 -Inf 的数量

data = subset(data, select = c("gene_id",pval_col_name))
data = data[order(data[,pval_col_name]), , drop = FALSE]
data = data[!duplicated(data[["gene_id"]]), , drop = FALSE]

# Function to compute q-values with fallback on error
compute_qvalues <- function(pval_col_name, data) {
  # Ensure the p-values are in the range [0, 1]
  # Initialize qobj
  qobj <- NULL
  
  # Try to compute q-values using bootstrap method
  tryCatch({
    qobj <- qvalue(data[[pval_col_name]])
  }, error = function(e) {
    message("Error in default qvalue method: ", e$message)
    # Fallback to simple pi0 estimation method if bootstrap fails
    message("Using simple pi0 estimation")
    pi0 <- sum(data[[pval_col_name]] > 0.05) / length(data[[pval_col_name]])
    qobj <<- qvalue(data[[pval_col_name]], pi0 = pi0)
  })
  
  if (is.null(qobj)) {
    stop("Failed to compute q-values using both methods.")
  }
  
  return(qobj)
}

qobj <- compute_qvalues(pval_col_name, data)

threshold_idx = which.min(abs(qobj$qvalues - fdr_threshold))
result <- data.frame (fdr  = qobj$qvalues[threshold_idx], pvalue = data[[pval_col_name]][threshold_idx])
print("FDR_result")
print(result)
write.table(result, output_file_path, sep = "\t", row.names = FALSE, col.names = TRUE)
