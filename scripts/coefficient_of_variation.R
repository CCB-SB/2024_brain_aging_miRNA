suppressMessages(library(data.table))
suppressMessages(library(openxlsx))


#snakemake = readRDS("snakemake_coefficient_of_variation.rds")

data.table::setDTthreads(1)
set.seed(snakemake@params$parameters_props$set_seed)

print("coefficient_of_variation")


#---------------------------------- Functions ----------------------------------
cv = function(mat){
  apply(mat, 1, function(x) sd(x, na.rm=T)) / rowMeans(mat, na.rm=T)
}


#---------------------------------- Load data ----------------------------------
data_input = snakemake@params$data_input
top_list = snakemake@params$top_list
results_folder = snakemake@params$results_folder

# save rdata
#saveRDS(snakemake, file = "snakemake_coefficient_of_variation.rds")

expr = fread(sprintf("%s_%s/%s_%s_quantification_%s.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_input$data_sub_set, data_input$rna_class, data_input$feature_filtering, data_input$norm, data_input$detection_rate, data_input$detection_group), sep='\t', header=T)
annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 

output_folder_table = sprintf("%s/%s_%s/results_%s/matrices/cv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
dir.create(output_folder_table, recursive=TRUE)


#------------------------------------ Script ----------------------------------- 
# save rnas
rna_ids = expr[[colnames(expr)[1]]]
# keep only expression
expr = expr[, sapply(expr, is.numeric), with=F]

# only expressed
expr_filtered = expr[rowSums(expr != 0) > 0,]
rna_ids_filtered = rna_ids[rowSums(expr != 0) > 0]

# calculate coefficient of variation
expr_filtered_cv = abs(cv(expr_filtered)) # for every row (miRNA) std/mean
ranked_filtered_cv = frank(-expr_filtered_cv, ties.method = "first") #-> Standardabweichung groß (und Mittelwert klein),  -> Standabweichung klein (und Mittelwert groß)

for(i in top_list){
  row_names = rna_ids_filtered[ranked_filtered_cv <= i]
  row_names_df = data.frame(feature = row_names)
  
  colnames(row_names_df) = data_input$rna_class
  
  # save pvca values for all tissues 
  fwrite(row_names_df, sprintf("%s/%s.csv", output_folder_table, i), sep = "\t", row.names = FALSE)
  write.xlsx(row_names_df, sprintf("%s/%s.xlsx", output_folder_table, i), colNames = TRUE, rowNames = FALSE, append = FALSE)
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
