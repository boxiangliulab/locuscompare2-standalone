# input file first column should be "rank", rest columns should be gene ranking of different tools, the order of the genes is used as the ranking
# method are one of: GEO, stuart
# results are p-value, the smaller the better
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("At least 2 arguments are rquired: input_file_path, output_file_path")
}
input_file_path = args[1]
output_file_path = args[2]

if (!is.na(args[3])) {
  sample_zie = as.integer(args[3])
} else {
  sample_zie = NA
}

method = args[4]
if (is.na(method)  || tolower(method) == 'na' || tolower(method) == 'none' || tolower(method) == "geo") {
  method = "geom.mean"
}

if (is.na(input_file_path) || !file.exists(input_file_path) || !file.info(input_file_path)$size > 0) {
  print("input_file_path file does not exist or is empty!")
  q(save = "no")
}

if (!require(RobustRankAggreg)) {
  stop("RobustRankAggreg not installed")
}

# row.names can not have duplicates. no other columns are allowed
data = read.table(file = input_file_path, row.names = "rank", sep = "\t", header = TRUE)
glist = list()
for (i in 1:ncol(data)) {
  col = data[i]
  # NA is not allowd in individual rank list; rank list may not be the same length
  col = col[!is.na(col)]
  glist = append(glist, list(col))
}
rra_result = aggregateRanks(glist = glist, method = method, N = sample_zie)
colnames(rra_result) = c("gene_id", paste(if (method == "geom.mean") "geo" else method, "ranking", sep = "_"))
# write result to output file
write.table(rra_result, if (endsWith(output_file_path, ".gz")) gzfile(output_file_path) else output_file_path, sep = "\t", row.names = FALSE, col.names = TRUE)
