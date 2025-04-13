library(data.table)


#---------------------------------- Load data ----------------------------------
data_input = snakemake@params$data_input
min_detection_list = snakemake@params$min_detection_list


#------------------------------------ Script ----------------------------------- 
for (rna_class in data_input$rna_class_list) {
  
  for (min_detection in min_detection_list) {
    # load expr
    # get the second part of the data sub_set
    if (data_input$data_sub_set %in% c("microglia", "ROSMAP")) {
      suffix = data_input$data_sub_set
    } else {
      parts = strsplit(data_input$data_sub_set, "_")[[1]]
      suffix = paste(parts[-1], collapse = "_")
    }
    tbl = fread(sprintf("%s_%s/%s_expression_raw_%s.tsv", data_input$preprocesing_data_folder, data_input$data_sub_set, rna_class, suffix), sep='\t', header=T)
    #print(tbl)
    sample_cols = colnames(tbl)[2:ncol(tbl)]
    tbl[, (sample_cols):=lapply(.SD, function(x) as.numeric(x >= min_detection)), .SDcols = sample_cols]

    #print(sprintf("%s%s/%s_complete_detection_matrix_min_detect=%s.csv", data_input$preprocesing_data_folder, data_input$data_sub_set, rna_class, min_detection))
    fwrite(tbl, sprintf("%s_%s/%s_complete_detection_matrix_min_detect=%s.csv", data_input$preprocesing_data_folder, data_input$data_sub_set, rna_class, min_detection), sep='\t', na="NA")
  }
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
