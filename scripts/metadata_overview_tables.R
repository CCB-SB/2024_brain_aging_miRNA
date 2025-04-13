library(data.table)


#snakemake = readRDS("snakemake_metadata_overview_tables.rds")

print("sig_tissue_specific_miRNAs_line_plot")


#---------------------------------- Functions ----------------------------------
create_pooling = function(row_info, row_list, pooling) {

  names(pooling) = row_list
  annot[[sprintf("%s_pooling", row_info)]] = annot[[row_info]]
  for (elem in row_list) {
    annot[[sprintf("%s_pooling", row_info)]] = gsub(elem, pooling[[paste(elem)]], annot[[sprintf("%s_pooling", row_info)]])
  }
  
  return(annot)  
}


#---------------------------------- Load data ----------------------------------
data_input = snakemake@params$data_input
combinations = snakemake@params$combinations
results_folder = snakemake@params$results_folder

# save rdata
#saveRDS(snakemake, file = "snakemake_metadata_overview_tables.rds")

annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 

output_folder = sprintf("%s_%s/sample_overview", data_input$raw_data_folder, data_input$data_sub_set)
dir.create(output_folder, recursive=TRUE)


#------------------------------------ Script ----------------------------------- 
for(i in 1:length(combinations)){
  
  count_variable = names(combinations)[i]
  row_info = combinations[[i]]$g[1]
  col_info = combinations[[i]]$g[2]
  pooling = combinations[[i]]$pooling
  
  if ((col_info == "time") | (col_info == "timepoint") | (col_info == "age") | (col_info == "Time") | (col_info == "Timepoint") | (col_info == "Age")) {
    print(sprintf("%s get sorted", col_info))
    col_list = sort(unique(annot[[col_info]]))
  } else {
    col_list = unique(annot[[col_info]])
  }
  
  if ((row_info == "time") | (row_info == "timepoint") | (row_info == "age") | (row_info == "Time") | (row_info == "Timepoint") | (row_info == "Age")) {
    print(sprintf("%s get sorted", row_info))
    row_list = sort(unique(annot[[row_info]]))
  } else {
    row_list = unique(annot[[row_info]])
  }
  
  if (is.null(pooling) == FALSE) {
    print("pooling")
    print(row_list)
    print(pooling)
    annot = create_pooling(row_info, row_list, pooling)
    row_info = sprintf("%s_pooling", row_info)
    row_list = unique(annot[[row_info]])
  }
  
  info_table = c()
  for (row in row_list) {
    info_row_col = c()
    tmp = (annot[[row_info]] == row)
    annot_filtered = annot[tmp,]
    info_row_col[["all"]] = length(unique(annot_filtered[[count_variable]]))
    for (col in col_list) {
      tmp = (annot_filtered[[col_info]] == col)
      info_row_col[[paste(col)]] = length(unique(annot_filtered[tmp,][[count_variable]]))
    }
    info_table = rbind(info_table, info_row_col)
  }
  
  total_count_variable_elements = length(unique(annot[[count_variable]]))
  
  # overview_brain_region_age
  info_table_df = as.data.frame(info_table)
  rownames(info_table_df) = paste(row_list)

  fwrite(info_table_df, sprintf("%s/overview_%s_%s_%s_total_%s_%ss.csv", output_folder, row_info, col_info, count_variable, total_count_variable_elements, count_variable), sep = "\t", row.names = TRUE)
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
