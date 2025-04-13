suppressMessages(library(data.table))
suppressMessages(library(corrplot))
suppressMessages(library(viridisLite))
suppressMessages(library(pheatmap))
suppressMessages(library(DESeq2))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
suppressMessages(library(openxlsx))


data_input = snakemake@params$data_input
adjustment = snakemake@params$adjustment
methods = snakemake@params$method
results_folder = snakemake@params$results_folder

output_folder = sprintf("%s/%s_%s/results_%s/matrices/correlation_mirna_mrna", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
dir.create(output_folder, recursive=TRUE)

# load mrna expr
tbl = read.csv(sprintf("%s/data_%s/mRNA_normalized_counts.csv", data_input$raw_data_folder, data_input$data_sub_set), sep='\t', header=TRUE, check.names = FALSE)
names_expressed_mrna = rownames(tbl)
rownames(tbl) = c()
names_mrna = toupper(names_expressed_mrna)
expr_mrna = as.data.table(tbl)
print(sprintf("Shape of gene expression matrix: %s", paste(dim(expr_mrna), collapse = "x")))

# match columns
tbl = fread(sprintf("%%s/data_%s/quantification/%s_%s_quantification_%s.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_input$data_sub_set, data_input$rna_class, data_input$feature_filtering, data_input$norm, data_input$detection_rate, data_input$detection_group), sep='\t', header=T)
expr = tbl[, sapply(tbl, is.numeric), with=F]
names_mirna = tbl[[data_input$rna_class]]

# get intersection of column names
common_samples = intersect(colnames(expr_mrna), colnames(expr))

# select columns from intersection
expr_mrna = expr_mrna[, ..common_samples]
expr = expr[, ..common_samples]
print(dim(expr_mrna))
print(dim(expr))


# Set up parallel computing
n_cores = 64
# print(sprintf("running in %s threads", n_cores))
cl <- makeCluster(n_cores, type="FORK")
registerDoParallel(cl)

for (method in methods) {
  correlation = matrix(NA, nrow(expr), nrow(expr_mrna))
  pval = matrix(NA, nrow(expr), nrow(expr_mrna))
  
  for (i in 1:nrow(expr)) {
    x = as.numeric(expr[i, ])
    correlation[i, ] = foreach (j=1:nrow(expr_mrna), .combine = "rbind") %dopar% {
      # calculate correlation
      y = as.numeric(expr_mrna[j, ])
      cor(x, y, method=method)
    }
    pval[i, ] = foreach (j=1:nrow(expr_mrna), .combine = "rbind") %dopar% {
      # calculate correlation
      y = as.numeric(expr_mrna[j, ])
      cor.test(x, y, method=method)$p.value
    }
    if ((i %% 10) == 0) {
      print(i)
    }
  }
  
  stopCluster(cl)
  
  # adjustment
  padj = matrix(p.adjust(pval, adjustment), nrow(pval), ncol(pval))
  
  # correlation
  # add colnames and mirna column
  correlation = data.frame(correlation)
  colnames(correlation) = names_mrna
  correlation$rna = names_mirna
  # reverse column order such that mirna are in front
  correlation = correlation[, rev(colnames(correlation))]
  
  # pvalue
  pval = data.frame(pval)
  colnames(pval) = names_mrna
  pval$rna = names_mirna
  # reverse column order such that mirna are in front
  pval = pval[, rev(colnames(pval))]
  
  # adjusted
  padj = data.frame(padj)
  colnames(padj) = names_mrna
  padj$rna = names_mirna
  padj = padj[, rev(colnames(padj))]
  
  
  write.table(correlation, sprintf("%s/correlation_%s.csv", output_folder, method), sep='\t', row.names = F)
  write.table(pval, sprintf("%s/pval_%s.csv", output_folder, method), sep='\t', row.names = F)
  write.table(padj, sprintf("%s/padj_%s.csv", output_folder, method), sep='\t', row.names = F)
  
  write.xlsx(correlation, sprintf("%s/correlation_%s.xlsx", output_folder, method), colNames = TRUE, rowNames = FALSE, append = FALSE)
  write.xlsx(pval, sprintf("%s/pval_%s.xlsx", output_folder, method), colNames = TRUE, rowNames = FALSE, append = FALSE)
  write.xlsx(padj, sprintf("%s/padj_%s.xlsx", output_folder, method), colNames = TRUE, rowNames = FALSE, append = FALSE)
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
