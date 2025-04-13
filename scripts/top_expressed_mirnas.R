suppressMessages(library(data.table))
suppressMessages(library(openxlsx))
suppressMessages(library(matrixStats))

#snakemake = readRDS("snakemake_top_expressed_mirnas.rds")

data.table::setDTthreads(1)
set.seed(snakemake@params$parameters_props$set_seed)

print("top_expressed_mirnas")


#---------------------------------- Load data ----------------------------------
data_input = snakemake@params$data_input
top_list = snakemake@params$top_list
results_folder = snakemake@params$results_folder

# save rdata
#saveRDS(snakemake, file = "snakemake_top_expressed_mirnas.rds")

expr = fread(sprintf("%s_%s/%s_%s_quantification_%s.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_input$data_sub_set, data_input$rna_class, data_input$feature_filtering, data_input$norm, data_input$detection_rate, data_input$detection_group), sep='\t', header=T)

output_folder_table = sprintf("%s/%s_%s/results_%s/matrices/most_expressed", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
dir.create(output_folder_table, recursive=TRUE)


#------------------------------------ Script ----------------------------------- 
# save rnas
rna_ids = expr[[colnames(expr)[1]]]
# keep only expression
expr = expr[, sapply(expr, is.numeric), with=F]

# only expressed
expr_filtered = expr[rowSums(expr != 0) > 0,]
rna_ids_filtered = rna_ids[rowSums(expr != 0) > 0]

# calculate most expressed index
# calculate medians over all samples
medians_df = data.frame(feature=rna_ids_filtered, medians=rowMedians(as.matrix(expr_filtered)))

ranked_median_df = order(as.numeric(medians_df$medians), decreasing = TRUE)
ordered_medians_df = medians_df[ranked_median_df,]

for(i in top_list){
  if (i != "all") {
    i = as.numeric(i)
    # ordered from small to large
    # get the first i columns (rnas)
    top_medians_df = ordered_medians_df[1:i,]
    row_names_df = data.frame(feature = top_medians_df$feature)
    
    #print(max(tsi_df))
    #print(min(tsi_df))
    #print(top_tsi_df)
    
    # save names of the top features
    
    colnames(row_names_df)[1] = data_input$rna_class
    
    fwrite(row_names_df, sprintf("%s/top_%s.csv", output_folder_table, i), sep = "\t", row.names = FALSE)
    write.xlsx(row_names_df, sprintf("%s/top_%s.xlsx", output_folder_table, i), colNames = TRUE, rowNames = FALSE, append = FALSE)
  } else {
    # save names of all features
    
    tmp = data.frame(feature = rna_ids)
    colnames(tmp)[1] = data_input$rna_class
    fwrite(tmp, sprintf("%s/all.csv", output_folder_table), sep = "\t", row.names = FALSE)
    write.xlsx(tmp, sprintf("%s/all.xlsx", output_folder_table), colNames = TRUE, rowNames = FALSE, append = FALSE)
  }
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
