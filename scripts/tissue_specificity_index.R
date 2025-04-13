suppressMessages(library(data.table))
suppressMessages(library(openxlsx))
suppressMessages(library(matrixStats))

#snakemake = readRDS("snakemake_tissue_specificity_index.rds")

data.table::setDTthreads(1)
set.seed(snakemake@params$parameters_props$set_seed)

print("tissue_specificity_index")


#---------------------------------- Functions ----------------------------------
calculate_tsi = function(mat) {
  N = length(mat[,-1])
  sum_ti = 0
  for (ti in colnames(mat[,-1])) {
    sum_ti = sum_ti + ( 1 - (mat[[ti]] / max(mat[,-1])) )
  }
  return(sum_ti / (N - 1))
}


#---------------------------------- Load data ----------------------------------
data_input = snakemake@params$data_input
tissues = snakemake@params$tissues
top_list = snakemake@params$top_list
results_folder = snakemake@params$results_folder

# save rdata
#saveRDS(snakemake, file = "snakemake_tissue_specificity_index.rds")

expr = fread(sprintf("%s_%s/%s_%s_quantification_%s.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_input$data_sub_set, data_input$rna_class, data_input$feature_filtering, data_input$norm, data_input$detection_rate, data_input$detection_group), sep='\t', header=T)
annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 

output_folder_table = sprintf("%s/%s_%s/results_%s/matrices/tsi", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
dir.create(output_folder_table, recursive=TRUE)


#------------------------------------ Script ----------------------------------- 
# save rnas
rna_ids = expr[[colnames(expr)[1]]]
# keep only expression
expr = expr[, sapply(expr, is.numeric), with=F]

# only expressed
expr_filtered = expr[rowSums(expr != 0) > 0,]
rna_ids_filtered = rna_ids[rowSums(expr != 0) > 0]

# calculate tissue specificity index
# calculate medians over all organs
medians = data.frame(feature=rna_ids_filtered)
for (ti in unique(annot[[tissues]])) {
  ids = annot[annot[[tissues]] == ti,][[data_input$identifier_column]]
  tmp = colnames(expr) %in% ids
  medians[[ti]] = rowMedians(as.matrix(expr[, ..tmp]))
}

tsi = c()
for (mir in medians$feature) {
  medians_feature = medians[medians$feature == mir,]
  tsi[[mir]] = calculate_tsi(medians_feature)
}

tsi_df = data.frame(tsi, check.names = FALSE)

ranked_tsi_df = order(as.numeric(tsi_df), decreasing = TRUE)
ordered_tsi_df = tsi_df[, ranked_tsi_df]
# expr_filtered_cv = abs(cv(expr_filtered)) # for every row (miRNA) std/mean
# ranked_filtered_cv = frank(-expr_filtered_cv, ties.method = "first") #-> Standardabweichung groß (und Mittelwert klein),  -> Standabweichung klein (und Mittelwert groß)

for(i in top_list){
  # row_names = rna_ids_filtered[ranked_tsi_df <= i]
  # ordered from small to large
  # get the first i columns (rnas)
  top_tsi_df = ordered_tsi_df[, 1:i]
  row_names = colnames(top_tsi_df)
  row_names_df = data.frame(feature = row_names)
  colnames(row_names_df) = data_input$rna_class
  
  #print(max(tsi_df))
  #print(min(tsi_df))
  #print(top_tsi_df)
  
  # save pvca values for all tissues 
  fwrite(row_names_df, sprintf("%s/%s.csv", output_folder_table, i), sep = "\t", row.names = FALSE)
  write.xlsx(row_names_df, sprintf("%s/%s.xlsx", output_folder_table, i), colNames = TRUE, rowNames = FALSE, append = FALSE)
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
