suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(openxlsx))
#library(corrplot)
suppressPackageStartupMessages(library(Hmisc))

#snakemake = readRDS("snakemake_correlation_age_table.rds")

set.seed(snakemake@params$parameters_props$set_seed)

print("correlation_age_table")


#---------------------------------- Functions ----------------------------------


#---------------------------------- Load data ---------------------------------- 
data_input = snakemake@params$data_input
method = snakemake@params$method
adjustment = snakemake@params$adjustment
prop = snakemake@params$prop
results_folder = snakemake@params$results_folder
ID = data_input$identifier_column

# save rdata
#saveRDS(snakemake, file = "snakemake_correlation_age_table.rds")
#stop()

expr = fread(sprintf("%%s/data_%s/%s_%s_quantification_%s.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_input$data_sub_set, data_input$rna_class, data_input$feature_filtering, data_input$norm, data_input$detection_rate, data_input$detection_group), sep='\t', header=T)
annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 
fixed_cols = c(data_input$rna_class, annot[[data_input$identifier_column]])
expr = expr[, ..fixed_cols]

#------------------------------------ Script ----------------------------------- 
for (m in method) {
  time = prop
  
  output_folder_tab = sprintf("%s/%s_%s/results_%s/matrices/corr_plots/age", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
  dir.create(output_folder_tab, recursive=TRUE)
  
  correlation = c()
  corr_pvalues = c()

  timepoints = annot[[time]]
  if (!is.numeric(timepoints)) {
    unique_strings = unique(timepoints)
    string_to_number = setNames(seq_along(unique_strings), unique_strings)
    
    numeric_timepoints = string_to_number[timepoints]
  } else {
    numeric_timepoints = timepoints
  }

  correlation = as.data.frame(cor(t(expr[,2:dim(expr)[2]]), numeric_timepoints, method = m))
  # if all expr for every timepoint for a sample is 0 for a feature
  # then the correlation is NA
  # we then force the correlation to be 0
  na_cor = is.na(correlation)
  # uncomment this to show the corresponding lines in the expr_matrix
  # print(expr[,2:dim(expr)[2]][na_cor, ])
  # replace NA with 0
  correlation[na_cor, ] = 0
  # get pvalues from rcorr
  # $P gets p-values
  pvalues_tmp = as.data.frame(rcorr(t(expr[,2:dim(expr)[2]]), numeric_timepoints, type=m)$P)
  # select y because x=matrix, y=timepoints
  # and remove last value because it is correlation of timepoints with timepoints = NA
  corr_pvalues = head(pvalues_tmp$y, -1)
  # adjustment
  corr_pvalues = as.data.frame(p.adjust(corr_pvalues, method=adjustment))
  # remove all pvalues == NA and replace them with 1 (then it is not significant)
  na_pvalue = is.na(corr_pvalues)
  corr_pvalues[na_pvalue] = 1
  # check if anything is still NA
  #print(sprintf("NA in correlation: %s", any(is.na(correlation))))
  #print(sprintf("NA in pvalues: %s", any(is.na(corr_pvalues))))
  
  correlation_df = as.data.frame(correlation)
  corr_pvalues_df = as.data.frame(corr_pvalues)

  # save tables
  correlation_df = data.frame(feature = expr[[data_input$feature_column]], corr = unname(correlation_df))
  colnames(correlation_df) = c(gsub("feature", paste(data_input$rna_class), colnames(correlation_df)))
  fwrite(correlation_df, sprintf("%s/corr_method=%s_%ss_with_%s.csv", output_folder_tab, m, data_input$rna_class, time), sep = "\t", row.names = TRUE)
  write.xlsx(correlation_df, sprintf("%s/corr_method=%s_%ss_with_%s.xlsx", output_folder_tab, m, data_input$rna_class, time), colNames = TRUE, rowNames = TRUE, append = FALSE)
  
  corr_pvalues_df = data.frame(feature = expr[[data_input$feature_column]], adj_p_value = unname(corr_pvalues_df))
  colnames(corr_pvalues_df) = gsub("feature", paste(data_input$rna_class), colnames(corr_pvalues_df) )
  fwrite(corr_pvalues_df, sprintf("%s/corr_method=%s_adj_pvalues_%ss_with_%s.csv", output_folder_tab, m, data_input$rna_class, time), sep = "\t", row.names = TRUE)
  write.xlsx(corr_pvalues_df, sprintf("%s/corr_method=%s_adj_pvalues_%ss_with_%s.xlsx", output_folder_tab, m, data_input$rna_class, time), colNames = TRUE, rowNames = TRUE, append = FALSE)
  
  print(sprintf("Done: %s", m))
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
