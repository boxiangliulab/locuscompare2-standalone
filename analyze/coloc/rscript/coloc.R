# Assumed GWAS min columns: ["var_id_", "chrom", "position", "beta", "varbeta", "pvalue", "se"]
# Assumed eQTL min columns: ["var_id_", "chrom", "position", "beta", "varbeta", "pvalue", "se", "gene_id", "maf"]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("At least 5 arguments are rquired: output_file_path, gwas_path, eqtl_path, gwas_sample_size, eqtl_sample_size")
}
output_file_path = args[1]
gwas_path = args[2]
eqtl_path = args[3]
gwas_sample_size = as.integer(args[4])
eqtl_sample_size = as.integer(args[5])
gwas_type = args[6]
eqtl_type = args[7]
# Original source gwas file path to be appended to result file, optional
src_gwas_path = args[8]
# Original source eqtl file path to be appended to result file, optional
src_eqtl_path = args[9]
# prior probability a SNP is associated with trait 1, default 1e-4
tp1 = as.numeric(args[10])
# prior probability a SNP is associated with trait 2, default 1e-4
tp2 = as.numeric(args[11])
# prior probability a SNP is associated with both traits, default 1e-5
tp12 = as.numeric(args[12])
# The threshold to consider overall H4 is true, only when overall H4 is true, the SNP level H4 is considered as relevant, optional
overall_h4_threshold = as.numeric(args[13])

if (is.na(gwas_path) || !file.exists(gwas_path) || !file.info(gwas_path)$size > 0) {
  print("GWAS file does not exist or is empty!")
  q(save = "no")
}
if (is.na(eqtl_path) || !file.exists(eqtl_path) || !file.info(eqtl_path)$size > 0) {
  print("eQTL file does not exist or is empty!")
  q(save = "no")
}
if (is.na(gwas_sample_size)) {
  stop("must supply argument gwas_sample_size")
}
if (is.na(eqtl_sample_size)) {
  stop("must supply argument eqtl_sample_size")
}
gwas_df = read.table(file = gwas_path, header = T)
eqtl_df = read.table(file = eqtl_path, header = T)
if (is.na(gwas_type) || tolower(gwas_type) == 'na' || tolower(gwas_type) == 'none') {
  gwas_type = "cc"
}
if (is.na(eqtl_type) || tolower(eqtl_type) == 'na' || tolower(eqtl_type) == 'none') {
  eqtl_type = "quant"
}
if (is.na(overall_h4_threshold)) {
  overall_h4_threshold = 0
}

if (!require(coloc)) {
  stop("coloc not installed")
}
# Sort eqtl_df according gwas_df
eqtl_df = eqtl_df[match(gwas_df$var_id_, eqtl_df$var_id_),]
# snplist = scan(file = ld_snp_list_path, character(), sep = "\n")
# ld_matrix = data.matrix(read.table(file = ld_file_path, row.names = snplist, col.names = snplist))
# if (nrow(ld_matrix) != ncol(ld_matrix)) {
#   stop("LD matrix is not square")
# }
# # Remove redundant data from ld_matrix
# ld_matrix = ld_matrix[gwas_df$var_id_, gwas_df$var_id_]
# ld_matrix[is.nan(ld_matrix)] = 0
input = merge(gwas_df, eqtl_df, by = "var_id_", all = FALSE, suffixes = c("_gwas", "_eqtl"))
# if (length(setdiff(input$var_id_, colnames(ld_matrix)))) {
#   stop("colnames in LD matrix do not contain all common SNPs")
# }
if (nrow(input) == 0) {
  print("No common snps found in two input dataframes")
  q(save = "no")
}

d1 = list(snp = input$var_id_, beta = input$beta_gwas, varbeta = input$varbeta_gwas, position = input$position_gwas, type = gwas_type, N = gwas_sample_size, MAF = input$maf)
d2 = list(snp = input$var_id_, beta = input$beta_eqtl, varbeta = input$varbeta_eqtl, position = input$position_eqtl, type = eqtl_type, N = eqtl_sample_size, MAF = input$maf)
# runsusie without extra arguments may occur error: "The estimated prior variance is unreasonably large"
# s1 = runsusie(d1)
# s2 = runsusie(d2)

# Specify prior_variance, either a scalar, or a vector of length L
# s1 = runsusie(d1, prior_variance = d1_prior_variance, estimate_prior_variance = FALSE)
# s2 = runsusie(d2, prior_variance = d2_prior_variance, estimate_prior_variance = FALSE)

# Use default prior_variance
# s1 = try(runsusie(d1, estimate_prior_variance = FALSE, maxit = 10000, repeat_until_convergence = FALSE))
# if (class(s1) == "try-error") {
#   q(save = "no")
# }
# s2 = try(runsusie(d2, estimate_prior_variance = FALSE, maxit = 10000, repeat_until_convergence = FALSE))
# if (class(s2) == "try-error") {
#   q(save = "no")
# }

# check_prior = FALSE may cause susieR not converge and then keep calculating
# s1 = runsusie(d1, check_prior = FALSE)
# s2 = runsusie(d2, check_prior = FALSE)

# result = coloc.susie(s1, s2)

result = coloc.abf(dataset1 = d1, dataset2 = d2, p1 = tp1, p2 = tp2, p12 = tp12)

if (is.null(result$summary)) {
  print(paste("No results for gene", input$gene_id[0:1]))
  q(save = "no")
}

# filter rows by H4, and then write to file
# summary_set = result$summary[result$summary$PP.H4.abf > 0.95,]
summary_set = result$summary

if (summary_set["PP.H4.abf"] < overall_h4_threshold) {
  print(paste("overall H4 is", summary_set["PP.H4.abf"], "for gene", input$gene_id[0:1], "which is too small and will not be included in the result file"))
  q(save = "no")
}

# if (nrow(summary_set) == 0) {
#   print(paste("No records with H4 > specified cut off found for gene ", input$gene_id[0:1]))
#   q(save = "no")
# }
appendMode = file.exists(output_file_path) && file.info(output_file_path)$size > 0
# append extra info to result
extra_info_df = input[, c("var_id_", "chrom_gwas", "gene_id")]
colnames(extra_info_df)[which(colnames(extra_info_df) == "chrom_gwas")] = "chrom"
detail_result = result$results
detail_result$overall_H4 = summary_set["PP.H4.abf"]
if (!is.na(src_gwas_path)) {
  detail_result$gwas_path = src_gwas_path
}
if (!is.na(src_eqtl_path)) {
  detail_result$eqtl_path = src_eqtl_path
}
detail_result_with_gene = merge(detail_result, extra_info_df, by.x = "snp", by.y = "var_id_", all.x = TRUE)
write.table(detail_result_with_gene, if (endsWith(output_file_path, ".gz")) gzfile(output_file_path) else output_file_path, append = appendMode, sep = "\t", row.names = FALSE, col.names = !appendMode)
