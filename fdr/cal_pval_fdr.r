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
# 将pval_col_name这一列的NA值替换为1
data[[pval_col_name]][is.na(data[[pval_col_name]])] <- 1
# 确保pval_col_name列的值在[0, 1]范围内
data = data %>% filter(data[[pval_col_name]] >= 0 & data[[pval_col_name]] <= 1)

data = subset(data, select = c("gene_id",pval_col_name))
data = data[order(data[,pval_col_name]), , drop = FALSE]
data = data[!duplicated(data[["gene_id"]]), , drop = FALSE]
qobj = qvalue(data[[pval_col_name]])
threshold_idx = which.min(abs(qobj$qvalues - fdr_threshold))
result <- data.frame (fdr  = qobj$qvalues[threshold_idx], pvalue = data[[pval_col_name]][threshold_idx])
write.table(result, output_file_path, sep = "\t", row.names = FALSE, col.names = TRUE)
