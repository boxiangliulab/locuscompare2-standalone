# input file first column should be "gene_id", other required columns: probability, zscore, coloumns other than these are ignored
# results are probability, the larger the better
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("At least 2 arguments are rquired: input_file_path, output_file_path")
}
input_file_path = args[1]
output_file_path = args[2]
prior_fun = args[3]
if (is.na(prior_fun)) {
  prior_fun = "expit"
}
if (!file.exists(input_file_path) || !file.info(input_file_path)$size > 0) {
  print("input_file_path file does not exist or is empty!")
  q(save = "no")
}

if (!require(INTACT)) {
  stop("INTACT not installed")
}

# row.names can not have duplicates
data = read.table(file = input_file_path, row.names = "gene_id", sep = "\t", header = TRUE)
# run intact
data["intact_probability"] = intact(GLCP_vec = data[["probability"]], prior_fun = get(prior_fun), z_vec = data[["zscore"]], t = 0.01)
# add gene_id column and move it to the first column
data["gene_id"] = rownames(data)
cols = colnames(data)
col_order = c("gene_id", cols[!cols %in% c("gene_id")])
data = data[, col_order]
# sort dataframe by intact result
data = data[order(data["intact_probability"], decreasing = TRUE),]
# write result to output file
write.table(data, if (endsWith(output_file_path, ".gz")) gzfile(output_file_path) else output_file_path, sep = "\t", row.names = FALSE, col.names = TRUE)
