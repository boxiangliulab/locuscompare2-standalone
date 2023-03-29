# input file must contain "gene_id" column
# results are tab seperated table with 2 column (1 row): fdr, pvalue
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("At least 2 arguments are rquired: input_file_path, output_file_path")
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
data = subset(data, select = c("gene_id",pval_col_name))
data = data[order(data[pval_col_name]), , drop = FALSE]
data = data[!duplicated(data[["gene_id"]]), , drop = FALSE]
fdr = p.adjust(data[[pval_col_name]], "BH")
threshold_idx = which.min(abs(fdr - fdr_threshold))
result <- data.frame (fdr  = fdr[threshold_idx], pvalue = data[[pval_col_name]][threshold_idx])
write.table(result, output_file_path, sep = "\t", row.names = FALSE, col.names = TRUE)
