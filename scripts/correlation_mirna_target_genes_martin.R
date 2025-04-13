suppressMessages(library(data.table))
suppressMessages(library(DESeq2))
suppressMessages(library(openxlsx))
suppressMessages(library(foreach))

#snakemake = readRDS("snakemale_files/snakemake_correlation_mirna_target_genes_martin.rds")

set.seed(snakemake@params$parameters_props$set_seed)

print("correlation_mirna_target_genes_martin")


# ---------------------- Load data ----------------------
data_input = snakemake@params$data_input
feature = snakemake@params$feature
gene_list_cas_folder_path = snakemake@params$gene_list_cas_folder_path
gene_list_target_genes = snakemake@params$gene_list_target_genes
split_prop = snakemake@params$split_prop
corr_params = snakemake@params$corr_params
results_folder = snakemake@params$results_folder

# save rdata
#saveRDS(snakemake, file = "snakemake_correlation_mirna_target_genes_martin.rds")
#stop()

# load annot
annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character"))

# load mrna expr
tbl = read.csv(sprintf("%s_%s/mRNA_normalized_counts.csv", data_input$raw_data_folder, data_input$data_sub_set), sep='\t', header=TRUE, check.names = FALSE)
names_expressed_mrna = rownames(tbl)
rownames(tbl) = c()
names_mrna = toupper(names_expressed_mrna)
expr_mrna = as.data.table(tbl)
print(sprintf("Shape of gene expression matrix: %s", paste(dim(expr_mrna), collapse = "x")))

# load mirna expr
tbl = fread(sprintf("%s_%s/%s_%s_quantification_%s.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_input$data_sub_set, data_input$rna_class, data_input$feature_filtering, data_input$norm, data_input$detection_rate, data_input$detection_group), sep='\t', header=T)
expr_mirna = tbl[, sapply(tbl, is.numeric), with=F]
names_mirna = tbl[[data_input$rna_class]]
print(sprintf("Shape of miRNA expression matrix: %s", paste(dim(expr_mirna), collapse = "x")))

# load cas list
gene_list_cas = fread(gene_list_cas_folder_path, sep='\t', header=T)
gene_list_cas$Gene = toupper(gene_list_cas$Gene)
print(sprintf("Length of cas gene list: %s", length(gene_list_cas$Gene)))

# load target gene list
#gene_list_target_genes_all = fread(gene_list_target_genes_folder_path, sep='\t', header=T)
#gene_list_target_genes_all$`Target / Target pathway` = toupper(gene_list_target_genes_all$`Target / Target pathway`)
#print(sprintf("Length of unique target gene list: %s", length(unique(gene_list_target_genes_all$`Target / Target pathway`))))
#gene_list_target_genes_functional = gene_list_target_genes_all[gene_list_target_genes_all$Support == "Functional MTI",]
#print(sprintf("Length of unique target gene list without weak: %s", length(unique(gene_list_target_genes_functional$`Target / Target pathway`))))

# load supplementary table from paper 10.1080/15476286.2025.2449775 Hart et al.
tmp = read.xlsx(gene_list_target_genes$path, sheet = gene_list_target_genes$sheet, startRow = 2, colNames = TRUE, rowNames = FALSE)$Hart.et.al.
gene_list_target_genes_martin = as.vector(na.omit(tmp))
print(sprintf("Length of unique target gene from paper: %s", length(unique(gene_list_target_genes_martin))))

# create output folders
output_folder = sprintf("%s/%s_%s/results_%s/matrices/%s_target_genes_martin_corr", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, feature)
dir.create(output_folder, recursive=TRUE)


#------------------------------------ Script ----------------------------------- 
# take eather all target genes or only non weak ones
#i = 1
#for (gene_list_target_genes in list(unique(gene_list_target_genes_all$`Target / Target pathway`), unique(gene_list_target_genes_functional$`Target / Target pathway`))) {
#  if (i == 1) {
#    target_genes_source = "target_genes_all"
#  } else {
#    target_genes_source = "target_genes_without_weak"
#  }
  
target_genes_source = "target_genes_martin"

# is there an overlap between the cas genes and the target genes
common_genes = intersect(gene_list_cas$Gene, gene_list_target_genes_martin)
print(sprintf("Between CAS and the target genes from the paper, there is an overlap of %s genes (%s)", length(common_genes), paste(common_genes, collapse = ", ")))

# filter samples
# get intersection of column names
common_samples = intersect(colnames(expr_mrna), colnames(expr_mirna))

# select columns from intersection
expr_mrna = expr_mrna[, ..common_samples]
expr_mirna = expr_mirna[, ..common_samples]
print(sprintf("Shape after intersection filtering of gene expression matrix: %s", paste(dim(expr_mrna), collapse = "x")))
print(sprintf("Shape after intersection filtering of miRNA expression matrix: %s", paste(dim(expr_mirna), collapse = "x")))

# filter annot also to only keep the samples which are in both expression matrices
print(sprintf("Number of samples in annot before filtering: %s", dim(annot)[1]))
tmp = (annot[[data_input$identifier_column]] %in% common_samples)
annot_intersection = annot[tmp,]
print(sprintf("Shape of annot after filtering: %s", dim(annot_intersection)[1]))

# filter features
# filter genes by target gene list
tmp = (names_mrna %in% gene_list_target_genes_martin)
expr_mrna_genes = expr_mrna[tmp,]
names_mrna_genes = names_mrna[tmp]

# if there are some genes not included in the expression file
if (length(setdiff(gene_list_target_genes_martin, names_mrna)) > 0) {
  print(sprintf("There are %s (%s) genes missing in the expression which a contained in the inserted target gene file", length(setdiff(gene_list_target_genes_martin, names_mrna)), paste(setdiff(gene_list_target_genes_martin, names_mrna), collapse = ", ")))
}
print(sprintf("Shape of gene expression matrix: %s", paste(dim(expr_mrna_genes), collapse = "x")))

# filter miRNA expression by given feature from the config
expr_mirna_feature = expr_mirna[names_mirna == feature,]

correlation_ti = c()
padj_ti = c()
for (ti in unique(annot_intersection[[split_prop]])) {
  # filter expression files for tissue ti
  group_ids = annot_intersection[annot_intersection[[split_prop]] == ti,][[data_input$identifier_column]]
  expr_mrna_genes_ti = expr_mrna_genes[,..group_ids]
  expr_mirna_feature_ti = expr_mirna_feature[,..group_ids]
  if (all(colnames(expr_mrna_genes_ti) == colnames(expr_mirna_feature_ti))) {
  } else {
    print(ti)
  }
  
  # correlation calculation   
  for (method in corr_params$method) {
    #print(sprintf("%s: %s", method, ti))
    # create folder
    dir.create(sprintf("%s/%s/%s", output_folder, target_genes_source, method), recursive=TRUE)
    
    correlation = matrix(NA, nrow(expr_mirna_feature_ti), nrow(expr_mrna_genes_ti))
    pval = matrix(NA, nrow(expr_mirna_feature_ti), nrow(expr_mrna_genes_ti))
    
    for (i in 1:nrow(expr_mirna_feature_ti)) {
      x = as.numeric(expr_mirna_feature_ti[i, ])
      correlation[i, ] = foreach (j=1:nrow(expr_mrna_genes_ti), .combine = "rbind") %do% {
        # calculate correlation
        y = as.numeric(expr_mrna_genes_ti[j, ])
        cor(x, y, method=method)
      }
      pval[i, ] = foreach (j=1:nrow(expr_mrna_genes_ti), .combine = "rbind") %do% {
        # calculate correlation
        y = as.numeric(expr_mrna_genes_ti[j, ])
        cor.test(x, y, method=method)$p.value
      }
      if ((i %% 10) == 0) {
        print(i)
      }
    }
    
    # adjustment
    padj = matrix(p.adjust(pval, corr_params$adjustment), nrow(pval), ncol(pval))
    
    # correlation
    # add colnames and mirna column
    correlation = data.frame(correlation)
    colnames(correlation) = names_mrna_genes
    correlation$rna = feature
    # reverse column order such that mirna are in front
    correlation = correlation[, rev(colnames(correlation))]
    
    # pvalue
    pval = data.frame(pval)
    colnames(pval) = names_mrna_genes
    pval$rna = feature
    # reverse column order such that mirna are in front
    pval = pval[, rev(colnames(pval))]
    
    # adjusted
    padj = data.frame(padj)
    colnames(padj) = names_mrna_genes
    padj$rna = feature
    padj = padj[, rev(colnames(padj))]

    # save corr data files
    write.table(correlation, sprintf("%s/%s/%s/correlation_%s.csv", output_folder, target_genes_source, method, ti), sep='\t', row.names = F)
    write.table(pval, sprintf("%s/%s/%s/pval_%s.csv", output_folder, target_genes_source, method, ti), sep='\t', row.names = F)
    write.table(padj, sprintf("%s/%s/%s/padj_%s.csv", output_folder, target_genes_source, method, ti), sep='\t', row.names = F)
    
    write.xlsx(correlation, sprintf("%s/%s/%s/correlation_%s.xlsx", output_folder, target_genes_source, method, ti), colNames = TRUE, rowNames = FALSE, append = FALSE)
    write.xlsx(pval, sprintf("%s/%s/%s/pval_%s.xlsx", output_folder, target_genes_source, method, ti), colNames = TRUE, rowNames = FALSE, append = FALSE)
    write.xlsx(padj, sprintf("%s/%s/%s/padj_%s.xlsx", output_folder, target_genes_source, method, ti), colNames = TRUE, rowNames = FALSE, append = FALSE)
    
    correlation_ti[[ti]] = correlation
    padj_ti[[ti]] = padj
  }
}
#  i = i + 1
#}
for (method in corr_params$method) {
  for (corr_th in c(0.5, 0.3)) {
    sig_neg_correlated_genes = c()
    info_table = c()
    for (ti in unique(annot_intersection[[split_prop]])){
      # if padj is not sorted like correlation we do so
      padj = padj_ti[[ti]][colnames(correlation_ti[[ti]])]
      
      # filter for sig. neg. correlated genes
      keep_genes = ((correlation_ti[[ti]] <= -corr_th) & (padj < corr_params$sig_lvl))
      sig_neg_correlated_genes[[ti]] = colnames(correlation_ti[[ti]][,keep_genes, drop = FALSE])
      
      for (rna in sig_neg_correlated_genes[[ti]]) {
        info_table = rbind(info_table, c(rna, ti, correlation_ti[[ti]][[rna]], padj[[rna]]))
      }
    }
    info_table_df = data.frame(info_table)
    colnames(info_table_df) = c("mRNA", split_prop, "corr_value", "adj_p_value")
    
    #sig_neg_correlated_genes_rna = c()
    info_table_rna = c()
    if (corr_th == 0.5) {
      for (r in info_table_df$mRNA) {
        for (ti in unique(annot_intersection[[split_prop]])) {
          # if padj is not sorted like correlation we do so
          padj = padj_ti[[ti]][colnames(correlation_ti[[ti]])]
          
          # filter for sig. neg. correlated genes
          keep = ((correlation_ti[[ti]][r] <= -0.3) & (padj[r] < corr_params$sig_lvl))
          #sig_neg_correlated_genes_rna[[ti]] = colnames(correlation_ti[[ti]][,keep_genes, drop = FALSE])
          if (keep) {
            info_table_rna = rbind(info_table_rna, c(r, ti, correlation_ti[[ti]][[r]], padj[[r]]))
          }
        }
      }
      info_table_rna_df = data.frame(info_table_rna)
      colnames(info_table_rna_df) = c("mRNA", split_prop, "corr_value", "adj_p_value")
    
      write.table(info_table_rna_df, sprintf("%s/%s/%s/sig_neg_correlation_corr_the=%s_adj.csv", output_folder, target_genes_source, method, corr_th), sep='\t', row.names = F)
      write.xlsx(info_table_rna_df, sprintf("%s/%s/%s/sig_neg_correlation_corr_the=%s_adj.xlsx", output_folder, target_genes_source, method, corr_th), colNames = TRUE, rowNames = FALSE, append = FALSE)
    }
    write.table(info_table_df, sprintf("%s/%s/%s/sig_neg_correlation_corr_the=%s.csv", output_folder, target_genes_source, method, corr_th), sep='\t', row.names = F)
    write.xlsx(info_table_df, sprintf("%s/%s/%s/sig_neg_correlation_corr_the=%s.xlsx", output_folder, target_genes_source, method, corr_th), colNames = TRUE, rowNames = FALSE, append = FALSE)
  }
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
